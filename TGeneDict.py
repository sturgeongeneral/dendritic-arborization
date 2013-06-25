#infile = open("C:\Users\Becky\Desktop\worm\c_elegans.current_development.cds_transcripts.fa","r")
#
#line=infile.readline()
GeneNameDict = {}
#while line != "":
#	if line[0] == ">":
#		line = line.strip()
#		line = line.replace("gene=WBGene", "")
#		line = line.replace(">", "")
#		parts = line.split(" ")
#		GeneNameDict[parts[1].lower()] = parts[0]
#	line=infile.readline()
#	
#infile = open("C:\Users\Becky\Desktop\worm\c_elegans.current.cds_transcripts.fa","r")
#
#line=infile.readline()
#while line != "":
#	if line[0] == ">":
#		line = line.strip()
#		line = line.replace("gene=WBGene", "")
#		line = line.replace(">", "")
#		parts = line.split(" ")
#		GeneNameDict[parts[1].lower()] = parts[0]
#	line=infile.readline()
#	
infile = open("C:\Users\Becky\Desktop\worm\c_elegans.WS235.functional_descriptions.txt","r")

line=infile.readline()
while line != "":
	if line[0] == "=":
		line = infile.readline()
		line = line.strip()
		if not line:
			continue
		line = line.replace("WBGene", "")
		parts = line.split("\t")
		if "not known" not in parts[1]:
			GeneNameDict[parts[0].lower()] = parts[1]
	line=infile.readline()

#base = "tri-output\\prop-and-upregulated-not-dauer"
base = "prop-not-upregulated"
infile = open('C:\Users\Becky\Desktop\worm\\new-ven-diagram\output\\' + base,"r")
outfile = open('c:\Users\Becky\Desktop\worm\\new-ven-diagram\output\\' + base + '.genes.txt',"w")

line=infile.readline()
while line != "":
	line = line.strip()
	if line in GeneNameDict:
		#got a match
		
		outfile.write(GeneNameDict[line] + "\n") 
	else:
		#did not find a match
		outfile.write("What is" + line + "?" + "\n") 
	line=infile.readline()