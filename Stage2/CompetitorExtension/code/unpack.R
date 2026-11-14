# Run from Stage2/CompetitorExtension. Each archive unpacks into its own directory.
args<-commandArgs(TRUE)
section<-match.arg(args[1],c('section4','section5','section6'))
dir.create(file.path('replay',section),recursive=TRUE,showWarnings=FALSE)
unzip(file.path('archives',paste0(section,'.zip')),exdir=file.path('replay',section))
