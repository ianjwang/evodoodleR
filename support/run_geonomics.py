"""Run Geonomics; export data for R. No evodoodle or plotting calls."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import random
os.environ.setdefault('MPLBACKEND', 'Agg')
import numpy as np
import geonomics as gnx


def progress(folder, message):
    # Atomic replacement prevents the UI from reading a partially written line.
    temporary = folder / 'progress.tmp'
    temporary.write_text(message + '\n')
    temporary.replace(folder / 'progress.txt')
    print(message, flush=True)


def export_snapshot(model, folder, stage, elapsed):
    progress(folder, f'Exporting data at step {elapsed}…')
    species = model.comm[0]
    prefix = folder / f'stage_{stage:02d}'
    ids = list(species)
    if not ids:
        raise RuntimeError(f'Population extinct at step {elapsed}; change the landscape or parameters.')
    # Explicit ID order keeps genotypes and locations aligned.
    genotypes = np.asarray(species._get_genotypes(individs=ids))
    if genotypes.ndim != 3 or genotypes.shape[2] != 2:
        raise RuntimeError(f'Unexpected diploid genotype shape: {genotypes.shape}')
    heterozygosity = np.mean(genotypes[:, :, 0] != genotypes[:, :, 1], axis=1)
    dosage = np.sum(genotypes, axis=2).astype(np.uint8)
    with gzip.open(str(prefix) + '_dosage.csv.gz', 'wt') as f:
        np.savetxt(f, dosage, fmt='%d', delimiter=',')
    with open(str(prefix) + '_individuals.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id', 'x', 'y', 'phenotype', 'heterozygosity'])
        for ident, het in zip(ids, heterozygosity):
            ind = species[ident]
            writer.writerow([ident, ind.x, ind.y, ind.z[0], het])
    with open(str(prefix) + '_coordinates.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['id', 'x', 'y'])
        for ident in ids:
            writer.writerow([ident, species[ident].x, species[ident].y])
    # Release plotting arrays before Geonomics constructs its VCF output.
    del genotypes, dosage
    model.write_gendata(str(prefix) + '.vcf', spp=0, n=None,
                        include_fixed_sites=True)
    return {'stage': stage, 'steps': elapsed, 'individuals': len(ids),
            'mean_heterozygosity': float(heterozygosity.mean())}


def run(folder, stages, phi, loci, factor, movement, radius, seed):
    folder = Path(folder).resolve()
    with open(folder / 'run_settings.json', 'w') as f:
        json.dump(dict(stages=stages, phi=phi, loci=loci, capacity_factor=factor,
                       movement=movement, mating_radius=radius, seed=seed), f, indent=2)
    np.random.seed(seed)
    random.seed(seed)
    with open(Path(__file__).with_name('model_params.json')) as f:
        params = json.load(f)
    layers = {name: np.loadtxt(folder / (name + '.csv'), delimiter=',', ndmin=2)
              for name in ('capacity', 'connectivity', 'environment')}
    shape = layers['capacity'].shape
    if any(a.shape != shape or not np.isfinite(a).all() or np.any((a < 0) | (a > 1))
           for a in layers.values()):
        raise ValueError('Landscapes must share dimensions and have finite values from 0 to 1.')
    params['landscape']['main']['dim'] = (shape[1], shape[0])
    for incoming, target in [('capacity', 'carrying_capacity'),
                             ('connectivity', 'connectivity'), ('environment', 'environment')]:
        params['landscape']['layers'][target]['init']['defined']['rast'] = layers[incoming]
    species = params['comm']['species']['spp_0']
    species['init'].update(K_layer='carrying_capacity', K_factor=factor)
    species['movement']['move_surf']['layer'] = 'connectivity'
    species['movement']['movement_distance_distr_param1'] = movement
    species['mating']['mating_radius'] = radius
    species['gen_arch']['traits']['trait_0'].update(layer='environment', phi=phi, n_loci=loci)
    params['model']['T'] = max(sum(stages), params['model']['T'])
    params['model']['num'] = seed
    # Keep incidental Geonomics outputs in the run folder.
    os.chdir(folder)
    progress(folder, 'Building the Geonomics model…')
    model = gnx.make_model(parameters=gnx.make_params_dict(params, 'evodoodleR'), verbose=True)
    progress(folder, 'Burning in the population (up to 10,000 steps)…')
    model.walk(T=10000, mode='burn')  # Same burn-in budget as the lab's original init_mod.
    if not model.comm.burned:
        raise RuntimeError('Burn-in did not finish within 10,000 steps.')
    summaries = [export_snapshot(model, folder, 0, 0)]
    elapsed = 0
    for i, steps in enumerate(stages, 1):
        progress(folder, f'Running main steps {elapsed + 1}–{elapsed + steps}…')
        model.walk(T=steps, mode='main')
        elapsed += steps
        summaries.append(export_snapshot(model, folder, i, elapsed))
    with open(folder / 'stages.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)
    progress(folder, 'Simulation complete; preparing R plots.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('folder')
    parser.add_argument('stages')
    parser.add_argument('phi', type=float)
    parser.add_argument('loci', type=int)
    parser.add_argument('factor', type=float)
    parser.add_argument('movement', type=float)
    parser.add_argument('radius', type=float)
    parser.add_argument('seed', type=int)
    a = parser.parse_args()
    run(a.folder, [int(s) for s in a.stages.split(',')], a.phi, a.loci,
        a.factor, a.movement, a.radius, a.seed)
