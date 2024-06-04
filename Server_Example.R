## Average Expression

output$tab = renderDataTable({ # Create table and use tab object in UI to show output
  fg = AverageExpression(data, assay = "RNA", group.by = "gast.id") # seurat function
  fg = fg$RNA # will do expression by RNA so $RNA is important
  fg = cbind(rownames(fg), fg)  $ # cbind helps show the rownames
  fg # return object 
})

output$dept <- downloadHandler( # for download
  filename = function(){"thename.csv"}, 
  content = function(fname){
    write.csv(fg, fname)
  }
)



# Differential Expression

cds.markers = reactive({ # reactive creates a reactive object important 
  input$RUN # create run button and use isolate function so the function does not run immediately 
  isolate(cds.markers = FindMarkers(data, ident.1 = input$celltype,assay = "RNA", group.by = "gast.id", logfc.threshold = input$logfc, min.pct = input$pct))
  cds.markers = cbind(rownames(cds.markers), cds.markers)
  cds.markers
})
output$diffyy = renderDataTable(cds.markers()) # use brackets with object as this is how you access the reactive object

output$down23 <- downloadHandler(
  filename = function(){"thename.csv"}, 
  content = function(fname){
    write.csv(cds.markers(), fname)
  }
)
