
Venn_diagram_module_UI<-function(id)
{
  ns<-NS(id)
  tagList(
    DT::dataTableOutput(ns("de_venn")),
    downloadButton(ns('download_venny'), 'Download plot'),
    plotOutput(ns("venn_diag"))
   
  )
}

Venn_diagram_module<-function(input,output,session,DE_genes,
                                      combination,wgcna_output)
{
  combo<-combination()
  result<-DE_genes()
 
#######venn diagram#####
#display interactive table that summarizes DE genes identified for comparisons

output$de_venn <- DT::renderDataTable({
  
    rows<-length(combo())
    modules<-NULL
    WGCNA_matrix<-NULL
    res<-data.frame(matrix(NA, nrow = length(combo()), ncol = 3))
    
    if(!is.null(wgcna_output()))
    {
      if((length(wgcna_output()$modules())>0))
      {
        mod<-wgcna_output()$modules()
        modules<-as.data.frame(table(mod))
        colnames(modules)<-c("Var1","numbers")
        entry<-c(combo(), levels(modules$Var1))
        rows<-length(entry)
        res<-data.frame(matrix(NA, nrow = rows, ncol = 3))
        colnames(res)<-c('Up regulated','Down regulated','Both')
        rownames(res)<-lapply(1:rows, function(i) {
         entry[[i]]
        
        })
       for(i in 1:length(combo()))
       {
        res[i,1]<-nrow(as.data.frame(result()[[i]][1]))
        res[i,2]<-nrow(as.data.frame(result()[[i]][2]))
        res[i,3]<-nrow(as.data.frame(result()[[i]][3]))
       }
       for(i in length(combo())+1:nrow(modules))
       {
        res[i,1:2]<-0
        res[i,3]<-result()[[i]][3]
       }
      }
    }
    else{
      rownames(res)<-lapply(1:length(combo()), function(i) {
        combo()[[i]]
        
      })
      colnames(res)<-c('Up regulated','Down regulated','Both')
      for(i in 1:length(combo()))
      {
        res[i,1]<-nrow(as.data.frame(result()[[i]][1]))
        res[i,2]<-nrow(as.data.frame(result()[[i]][2]))
        res[i,3]<-nrow(as.data.frame(result()[[i]][3]))
      }
    }
    DT::datatable(res,class = 'cell-border stripe',
                  selection = list(target = 'cell'),
                  extensions = list('Scroller'=NULL,'Buttons'=NULL),
                  options = list(deferRender = TRUE,scrollX = TRUE,scrollY = 150,scroller = TRUE,dom = 'Bfrtip',
                                 buttons = list(list(extend = 'collection',buttons = c('csv', 'excel', 'pdf'),
                                                             text = 'Download table'))),
                  escape = FALSE
                  
    )
})

#Display the venn diagram when user selects the sets to get the overlap of the sets
observeEvent(input$de_venn_cell_clicked,{
  
  selected <- input$de_venn_cells_selected
  venn_list<-list()
  comp<-lapply(1:length(combo()), function(i) {
    combo()[[i]]
    
  })
  
  if(!is.null(wgcna_output()))
  {
    if((length(wgcna_output()$modules())>0))
    {
      mod<-wgcna_output()$modules()
      modules<-as.data.frame(table(mod))
      colnames(modules)<-c("Var1","numbers")
      entry<-c(combo(), levels(modules$Var1))
      rows<-length(entry)
      entry<-c(as.vector(combo()), as.vector(modules$Var1))
      comp<-lapply(1:rows, function(i) {
        entry[i]
    })
    }
  }
  reg<-c('up regulated genes','down regulated genes','All DE genes')
  pal<-c('pink','light green','light blue','light yellow')
  col<-c()
  circ<-c()
  venn_names<-list()
  if(nrow(selected) %in% 2:4)
  {
    for(i in 1:nrow(selected))
    {
      row<-selected[i,1]
      col<-selected[i,2]
      genes<-NULL
      if(row>length(combo()))
      {
        mod<-wgcna_output()$modules()
        modules<-as.data.frame(table(mod))
        colnames(modules)<-c("Var1","numbers")
        WGCNA_matrix<-wgcna_output()$WGCNA_matrix()
        idx_w<-which(mod==modules$Var1[row-length(combo())])
        genes<-colnames(WGCNA_matrix)[idx_w]
      }
      else
      {
        df<-as.data.frame(result()[[row]][col])
        genes<-df[,1]
        rownames(df)<-genes
      }
      venn_list[[length(venn_list)+1]]<-genes
      venn_names[[length(venn_names)+1]]<-paste(comp[row]," ",reg[col])
      
    }
    alpha<-c()
    i<-nrow(selected)
    if(i==2)
    {
      col<-c(pal[1],pal[2])
      circ<-c(pal[1],pal[2])
      alpha<-c(0.5,0.5)
    }
    else if(i==3)
    {
      col<-c(pal[1],pal[2],pal[3])#c('lightpink1','lightgreen','lightgoldenrod1')
      circ<-c(pal[1],pal[2],pal[3])
      alpha<-c(0.5,0.5,0.5)
    }
    else if(i==4)
    {
      col<-c(pal[1],pal[2],pal[3],pal[4])#c('lightpink1','lightgreen','lightgoldenrod1','lightskyblue1')
      circ<-c(pal[1],pal[2],pal[3],pal[4])
      alpha<-c(0.5,0.5,0.5,0.5)
    }
    
    output$venn_diag <- renderPlot({
      venn.plot <- venn.diagram(venn_list, NULL, fill=col,alpha=alpha,col=circ,
                                cex = 2, cat.fontface=4, category.names=venn_names)
      grid.draw(venn.plot)                     
      
    })
    #download venn diagram
    output$download_venny <- downloadHandler(
      filename = paste(" Venn diagram of DE genes: ",nrow(selected),"comparisons",'.pdf'),
      content = function(file) {
        #png(file,width=800, height=500)
        pdf(file)
        venn.plot <-venn.diagram(venn_list, NULL, fill=col,alpha=alpha,col=circ,
                                 cex = 2, cat.fontface=4, category.names=venn_names)
        grid.draw(venn.plot) 
        dev.off()
        
      }
      
    )
  }
})
}
