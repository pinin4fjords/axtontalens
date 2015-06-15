library(Biostrings)

#########################################################
#
# Collect and check the inputs
#
#########################################################

cl = commandArgs(trailingOnly=T)
target_sequence <- cl[1]
monomers_file <- cl[2]
fragment_defs_file <- cl[3]

if (is.na(target_sequence) || is.na(fragment_defs_file) || ! file.exists(fragment_defs_file) || is.na(monomers_file) || ! file.exists(monomers_file)){
  stop("Usage: \"Rscript makeFragments.R <target sequence> <valid monomers file> <valid fragment definitions file>\"")
}

#########################################################
#
# Read the monomers
#
#########################################################

sequences <- readDNAStringSet(monomers_file)

#########################################################
#
# Read the fragment definitions
#
#########################################################

fragment_defs <- read.delim(fragment_defs_file, sep="\t", header=F, stringsAsFactors=F)
lengths <- fragment_defs[,1]
fragment_defs <- lapply(fragment_defs[,2], function(x) lapply(unlist(strsplit(x, '-')), function(y) unlist(strsplit(y, ','))) )
names(fragment_defs) <- as.character(lengths)

#########################################################
#
# Check that the sequence length is valid given the 
# fragment definititions
#
#########################################################

chars <- unlist(strsplit(target_sequence, ''))

if (! length(chars) %in% names(fragment_defs)){
  stop(paste("You must supply a sequence of between", min(as.numeric(names(fragment_pairs))), "and",  max(as.numeric(names(fragment_pairs)))))
}

fragment_def <- fragment_defs[[as.character(length(chars))]]

#########################################################
#
# Assemble the fragments from the sequence and the 
# fragment definitions
#
#########################################################

monomers <- unlist(fragment_def)

# First select the right monomers....

monomers_used <- split(
  unlist(lapply(1:length(monomers), function(m) paste(monomers[m], chars[m], sep=''))),
  unlist(lapply(1:length(fragment_def), function(fd) rep(fd, length(fragment_def[[fd]])) ))
)

# ... then pick the right sequences (as DNAString objects) 

fragments <- DNAStringSet(lapply(monomers_used, function(mu){
  do.call(xscat, sequences[mu])
}))
names(fragments) <- unlist(lapply(monomers_used, function(mu){
  paste(mu, collapse='-')
}))

#########################################################
#
# Ouptut to file
#
#########################################################

filename <- paste(target_sequence,"fa", sep='.')
writeXStringSet(fragments, file=filename)
print(paste("####### Sequences output to file", filename))