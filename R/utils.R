# This is really a function to check for a valid method
get_test_method = function(x) {
	if (method %in% names(SUPPORTED_METHODS)) {
		return(SUPPORTED_METHODS[[method]])
	} else if (method %in% names(HIDDEN_METHODS)) {
		return(HIDDEN_METHODS[[method]])
	} else {
		stop(sprintf("Error: invalid enrichment test requested: %s, contact developer.",method))
	}
}

# Recodes # of peaks to be:
# { 0, if no peaks
# { 1, if number of peaks is >= threshold argument
recode_peaks = function(num_peaks,threshold=1) {
  as.numeric(num_peaks >= threshold);
}

# Quick function to change column names of data frame
# (by name, not by position.)
change_names = function(dframe,name_list) {
  for (n in names(name_list)) {
    names(dframe)[names(dframe)==n] = name_list[[n]];
  }
  dframe;
}

# Remove a file extension from a file name.
strip_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  s = s[1:(length(s)-1)];
  s = paste(s,collapse=".");
  return(s);
}

# Get a file's extension.
get_ext = function(file) {
  s = unlist(strsplit(file,".",fixed=T));
  return(tail(s,1));
}

# Checks to see if all elements of a list argument are possible.
check_arg = function(arg,possible_values,value=F) {
  if (!all(arg %in% possible_values)) {
    if (value) {
      return(arg[!arg %in% possible_values]);
    } else {
      return(F);
    }
  } else {
    return(T);
  }
}

genome_to_organism = function(genome) {
  code = substr(as.character(genome),1,2);
  if (code == 'mm') {
    org = 'mmu';
  } else if (code == 'hg') {
    org = 'hsa';
  } else if (code == 'rn') {
    org = 'rno';
  } else if (code == 'dm') {
  	org = 'dme';
  }
  else {
    org = NULL;
  }

  if (is.null(org)) {
    stop("Error: genome requested is not supported.");
  }

  org;
}
