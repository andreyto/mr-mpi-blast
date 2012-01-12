#!/usr/bin/env python

import random
 

allDbNames = []
for i in range(0, 10):
    #print "refseq_genomic.0" + str(i)
    allDbNames.append("refseq_genomic.0" + str(i))
for i in range(10, 34):
    allDbNames.append("refseq_genomic." + str(i))

for i in range(0, 10):
    allDbNames.append("wgs.0" + str(i))
for i in range(10, 56):
    allDbNames.append("wgs." + str(i))
        
for i in range(0, 9):
    allDbNames.append("nt.0" + str(i))
    
for i in range(0, 10):
    allDbNames.append("htgs.0" + str(i))

print allDbNames, len(allDbNames)
random.shuffle(allDbNames)
print allDbNames, len(allDbNames)

dbChunkNameFileName = "dbchunks.txt"
dbChunkNameFile = open(dbChunkNameFileName, "w")
for i in range(0, len(allDbNames)):
    dbChunkNameFile.write(allDbNames[i])
    dbChunkNameFile.write('\n')
       

### 6 * 16 + 13, total 7 chunks * 16G
#start = 1
#end = 16
#for n in range(0,6):
    #aliasFileName = "DB.0" + str(n) + ".nal"
    #aliasFile = open(aliasFileName, "w")
    #aliasFile.write("#\n")
    #aliasFile.write("TITLE " + aliasFileName + '\n')
    #aliasFile.write("#\n")

    #DbFileList = allDbNames[start-1:end]
    #print "size allDbNames[%d:%d] = %d" % (start-1, end, len(allDbNames[start-1:end]))
    #aliasFile.write("DBLIST ")
    #for i in range(start-1, end):
        #aliasFile.write(allDbNames[i])
        #aliasFile.write(' ')
    #aliasFile.write("\n#")
    #aliasFile.close()
    #start += 16
    #end += 16

### The rest 13 files
#aliasFileName = "DB.06.nal"
#aliasFile = open(aliasFileName, "w")
#aliasFile.write("#\n")
#aliasFile.write("TITLE " + aliasFileName + '\n')
#aliasFile.write("#\n")
#start = 97
#end = 109
#DbFileList = allDbNames[start-1:end]
#print "size allDbNames[%d:%d] = %d" % (start-1, end, len(allDbNames[start-1:end]))
#aliasFile.write("DBLIST ")
#for i in range(start-1, end):
    #aliasFile.write(allDbNames[i])
    #aliasFile.write(' ')
#aliasFile.write("\n#")
#aliasFile.close()
 
## EOF
