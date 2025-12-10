# obtains the approximate rho to be used for shapers based on sample params

##############################################################
#             PROCESS COMMAND LINE ARGUMENTS                 #
##############################################################
options(error=traceback)

args = commandArgs(trailingOnly = TRUE)
print("Arguments received:")
print(args)
if(length(args) < 7) {stop("not enough Arguments received")}

nkl0 = as.numeric(args[1])
nk1 = as.numeric(args[2])
nk0 = as.numeric(args[3])
nl1 = as.numeric(args[4])
nl0 = as.numeric(args[5])
outputLoc = args[6]
shaPRSscriptLoc= args[7]

# load exteral functions
source(shaPRSscriptLoc) 


# nkl0 number of controls overlapping between studies
# nk1 number of cases in study k
# nk0 number of controls in study k
# nl1 number of cases in study l
# nl0 number of controls in study l
print(paste0("calculating rho from numSharedControls: ", nkl0, " / numCases1: ", nk1 , " / numControls1: ", nk0 , " / numCases2: ", nl1, " / numControls2: ", nl0))
rho = shaPRS_rho(nkl0 = nkl0, nk1 = nk1, nk0= nk0, nl1= nl1, nl0= nl0)


write.table(rho, outputLoc, sep = "\t", row.names = F, col.names = F, quote = FALSE)
print(paste("written rho (",rho ,") to",outputLoc))


