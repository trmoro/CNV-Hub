import os

f = open("XCNV/bin/XCNV","r")
d = f.read()
f.close()

n = open("XCNV.R","w")
n.write(d.replace("-e #!/usr/bin/Rscript","#!/usr/bin/Rscript") )
n.close()

os.system("chmod +x XCNV.R")
