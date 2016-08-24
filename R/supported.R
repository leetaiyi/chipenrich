supported_locusdefs = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  ldefs = grep("locusdef",data_files,value=T);
  unique(str_replace(ldefs,"locusdef\\.(.+?)\\.",""));
}

supported_read_lengths = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  mappas = grep("mappa",data_files,value=T);
  sort(unique(as.numeric(str_replace(str_extract(mappas,"(\\d+)mer"),"mer",""))));
}

supported_genesets = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  geneset_files = grep("geneset",data_files,value=T);
  unlist(unique(Map(function(x) x[[2]],str_split(geneset_files,"\\."))));
}

supported_genomes = function() {
  piqr = data(package = "chipenrich.data");
  data_files = piqr$results[,3];

  as.character(na.omit(unique(str_match(data_files,"(locusdef|mappa)\\.(\\w+)\\.")[,3])));
}

supported_methods = function() {
  return(names(SUPPORTED_METHODS));
}
