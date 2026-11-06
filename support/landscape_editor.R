# Shiny draws all three layers; Save returns a named list of R matrices.
draw_landscapes <- function(d=10, initial=NULL, title="Draw landscapes",
                            file=NULL, redraw=interactive()) {
  if (!is.null(file) && file.exists(file) && !redraw) return(readRDS(file))
  if (!interactive())
    stop("Run this drawing chunk interactively and save its landscapes before knitting.")
  if (!requireNamespace("shiny", quietly=TRUE) ||
      !requireNamespace("jsonlite", quietly=TRUE))
    stop("Install shiny and jsonlite before drawing.")
  keys <- c("carrying_capacity","connectivity","environment")
  grids <- setNames(rep(list(matrix(1,d,d)),3),keys)
  if (!is.null(initial)) {
    if (!is.list(initial) || !all(keys %in% names(initial)))
      stop("initial must contain carrying_capacity, connectivity, and environment matrices.")
    grids <- lapply(initial[keys],as.matrix)
  }
  good <- vapply(grids,function(x)
    is.numeric(x) && length(dim(x))==2 && all(dim(x)>=2) &&
    all(is.finite(x)) && all(x>=0 & x<=1),logical(1))
  if (!all(good) || !identical(dim(grids[[1]]),dim(grids[[2]])) ||
      !identical(dim(grids[[1]]),dim(grids[[3]])))
    stop("The three grids must have matching dimensions and values from 0 to 1.")
  nr <- nrow(grids[[1]]); nc <- ncol(grids[[1]])
  palettes <- list(carrying_capacity=c("#f7fbff","#c6dbef","#6baed6","#2171b5","#08306b"),
                   connectivity=c("#f7fcf5","#c7e9c0","#74c476","#238b45","#00441b"),
                   environment=c("#542788","#998ec3","#f7f7f7","#f1a340","#b35806"))
  config <- lapply(keys,function(key) list(id=key,rows=nr,cols=nc,
    values=as.vector(t(grids[[key]])),palette=palettes[[key]]))
  javascript <- paste0("const editorConfigs=",jsonlite::toJSON(config,auto_unbox=TRUE),";",'
  (() => {
    const editors={};
    editorConfigs.forEach(cfg=>{
      const values=cfg.values.slice(), original=values.slice();
      const canvas=document.getElementById("canvas_"+cfg.id);
      const ctx=canvas.getContext("2d"), w=canvas.width/cfg.cols, h=canvas.height/cfg.rows;
      let down=false, previous=null;
      const brush=document.getElementById("brush_"+cfg.id);
      const buttons=document.querySelectorAll("#buttons_"+cfg.id+" button");
      buttons.forEach(button=>button.addEventListener("click",()=>{
        brush.value=button.dataset.value;
        buttons.forEach(other=>other.setAttribute("aria-pressed",String(other===button)));
      }));
      function color(v){
        const t=Math.max(0,Math.min(1,v))*4,i=Math.min(3,Math.floor(t)),f=t-i;
        const a=cfg.palette[i],b=cfg.palette[i+1];
        const rgb=[1,3,5].map(k=>Math.round(parseInt(a.slice(k,k+2),16)*(1-f)+parseInt(b.slice(k,k+2),16)*f));
        return "rgb("+rgb.join(",")+")";
      }
      function draw() {
        for(let r=0;r<cfg.rows;r++) for(let c=0;c<cfg.cols;c++) {
          ctx.fillStyle=color(values[r*cfg.cols+c]);
          ctx.fillRect(c*w,r*h,w,h);
          ctx.strokeStyle="#aaa";ctx.lineWidth=0.5;ctx.strokeRect(c*w,r*h,w,h);
        }
      }
      function paint(e) {
        const box=canvas.getBoundingClientRect();
        const c=Math.floor((e.clientX-box.left)/box.width*cfg.cols);
        const r=Math.floor((e.clientY-box.top)/box.height*cfg.rows);
        if(c<0||c>=cfg.cols||r<0||r>=cfg.rows){previous=null;return;}
        const start=previous||[r,c],steps=Math.max(Math.abs(r-start[0]),Math.abs(c-start[1]),1);
        const value=Number(document.getElementById("brush_"+cfg.id).value);
        for(let i=0;i<=steps;i++){
          const rr=Math.round(start[0]+(r-start[0])*i/steps);
          const cc=Math.round(start[1]+(c-start[1])*i/steps);
          values[rr*cfg.cols+cc]=value;
        }
        previous=[r,c];draw();
      }
      canvas.addEventListener("pointerdown",e=>{
        e.preventDefault();down=true;previous=null;canvas.setPointerCapture(e.pointerId);paint(e);
      });
      canvas.addEventListener("pointermove",e=>{if(down)paint(e);});
      ["pointerup","pointercancel"].forEach(event=>canvas.addEventListener(event,()=>{down=false;previous=null;}));
      document.getElementById("reset_"+cfg.id).addEventListener("click",()=>{
        values.splice(0,values.length,...original);draw();
      });
      editors[cfg.id]=values;draw();
    });
    document.getElementById("save_grid").addEventListener("click",()=>{
      if(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected())
        Shiny.setInputValue("saved_grid",{layers:editors,nonce:Date.now()},{priority:"event"});
    });
  })();')
  labels <- c(carrying_capacity="Carrying capacity",connectivity="Connectivity",
              environment="Environment")
  descriptions <- c(carrying_capacity="0 = no capacity; 1 = maximum relative capacity",
                    connectivity="0 = low connectivity; 1 = high connectivity",
                    environment="0 to 1 = local trait optimum")
  cards <- lapply(keys,function(key) shiny::column(4,
    shiny::h3(labels[[key]]), shiny::p(descriptions[[key]]),
    shiny::tags$canvas(id=paste0("canvas_",key),width=420,height=420,
      style="width:100%;border:1px solid #777;touch-action:none"),
    shiny::tags$input(id=paste0("brush_",key),type="hidden",value="0"),
    shiny::tags$div(id=paste0("buttons_",key),class="paint-buttons",
      role="group",`aria-label`=paste(labels[[key]],"paint value"),
      lapply(seq(0,1,0.25),function(value) shiny::tags$button(
        type="button",`data-value`=as.character(value),
        `aria-pressed`=if(value==0) "true" else "false",
        as.character(value)))),
    shiny::tags$button(id=paste0("reset_",key),type="button","Reset this layer")))
  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(
      ".paint-buttons{display:flex;gap:6px;margin:12px 0}.paint-buttons button{flex:1;padding:8px 4px;border:1px solid #777;border-radius:4px;background:white;color:#222}.paint-buttons button[aria-pressed=true]{background:#23527c;color:white;border:2px solid #163954;font-weight:bold}.paint-buttons button:focus-visible{outline:3px solid #e59500;outline-offset:2px}"))),
    shiny::h2(title),
    shiny::p("Paint any of the three grids. Save once when all layers are ready."),
    do.call(shiny::fluidRow,cards),
    shiny::tags$hr(),
    shiny::tags$button(id="save_grid",type="button",class="btn btn-primary",
                       "Save all and return to notebook"),
    shiny::actionButton("cancel_editor","Cancel"),
    shiny::p("After saving, this window disconnects normally. Return to the notebook."),
    shiny::tags$script(shiny::HTML(javascript))
  )
  server <- function(input,output,session) {
    shiny::observeEvent(input$saved_grid,{
      payload <- input$saved_grid$layers
      shiny::req(all(keys %in% names(payload)))
      answer <- setNames(lapply(keys,function(key) {
        values <- as.numeric(unlist(payload[[key]]))
        shiny::req(length(values)==nr*nc,all(is.finite(values)),all(values>=0 & values<=1))
        matrix(values,nrow=nr,ncol=nc,byrow=TRUE)
      }),keys)
      if (!is.null(file)) {
        dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE)
        saveRDS(answer,file)
      }
      shiny::stopApp(answer)
    })
    shiny::observeEvent(input$cancel_editor,shiny::stopApp(NULL))
  }
  answer <- shiny::runApp(shiny::shinyApp(ui,server),
                         launch.browser=getOption("shiny.launch.browser",TRUE))
  if (is.null(answer)) stop("Drawing cancelled; no landscapes returned.")
  answer
}
