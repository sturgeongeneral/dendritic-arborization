import re
import os
import sys
from os.path import isfile, join


class States:
	LookingFor_ViewDetailInfo = 1
	LookingFor_ImgLongTermDauer1 = 2
	LookingFor_ImgLongTermDauer2 = 3
	LookingFor_ImgLongTermDauer3 = 4
	LookingFor_ImgLongTermDauer4 = 5
	LookingFor_SevenTabsTR1 = 6
	LookingFor_SevenTabsTR2 = 7
	LookingFor_SevenTabsTR3 = 8
	LookingFor_SevenTabsTR4 = 9
	LookingFor_DauerNumber1 = 10
	LookingFor_DauerNumber2 = 11
	LookingFor_DauerNumber3 = 12
	LookingFor_DauerNumber4 = 14
	LookingFor_Identifier = 15

class Gene:
	DauerNumber1 = ""
	DauerNumber2 = ""
	DauerNumber3 = ""
	DauerNumber4 = ""
	CDNAClone = ""
	WBGene = ""
	WB_Gene_ID = ""
	EnsemblOrg = ""
	LastChanceID = ""
	def hasIdentifier(self):
		return self.CDNAClone != "" or self.WBGene != "" or self.WB_Gene_ID != "" or self.EnsemblOrg != "" or self.LastChanceID != ""
	def MeetsConservativeReqs(self):
		return self.DauerNumber1 == "3" and self.DauerNumber2 == "3" and self.DauerNumber3 == "3" and self.DauerNumber4 == "3"
	def prettyOutput(self):  
		global IDLookup
		if self.WBGene != "":
			return self.WBGene + "\n"
		else:
			if self.CDNAClone and self.CDNAClone in IDLookup:
				return IDLookup[self.CDNAClone] + " !1\n"
			elif self.WB_Gene_ID and self.WB_Gene_ID in IDLookup:
				return IDLookup[self.WB_Gene_ID] + " !2\n"
			elif self.EnsemblOrg and self.EnsemblOrg in IDLookup:
				return IDLookup[self.EnsemblOrg] + " !3\n"
			elif self.LastChanceID and self.LastChanceID in IDLookup:
				return IDLookup[self.LastChanceID] + " !4\n"
			elif self.CDNAClone and keyPrefixSearch(self.CDNAClone, IDLookup):
				return IDLookup[keyPrefixSearch(self.CDNAClone, IDLookup)] + " !5 " + self.CDNAClone + " " + keyPrefixSearch(self.CDNAClone, IDLookup) + "\n"
			elif self.WB_Gene_ID and keyPrefixSearch(self.WB_Gene_ID, IDLookup):
				return IDLookup[keyPrefixSearch(self.WB_Gene_ID, IDLookup)] + " !6 " + self.WB_Gene_ID + " " + keyPrefixSearch(self.WB_Gene_ID, IDLookup) + "\n"
			elif self.EnsemblOrg and keyPrefixSearch(self.EnsemblOrg, IDLookup):
				return IDLookup[keyPrefixSearch(self.EnsemblOrg, IDLookup)] + " !7 " + self.EnsemblOrg + " " + keyPrefixSearch(self.EnsemblOrg, IDLookup) + "\n"
			elif self.LastChanceID and keyPrefixSearch(self.LastChanceID, IDLookup):
				return IDLookup[keyPrefixSearch(self.LastChanceID, IDLookup)] + " !8 " + self.LastChanceID + " " + keyPrefixSearch(self.LastChanceID, IDLookup) + "\n"
			else:
				return "?" + " WB_Gene_ID: " + self.WB_Gene_ID + \
					" CDNAClone: "+ self.CDNAClone + " Ensemble.org: " + self.EnsemblOrg + \
					" LastChanceID: " + self.LastChanceID + "\n"
		#return "WBGene: " + self.WBGene + " WB_Gene_ID: " + self.WB_Gene_ID + \
		#	" CDNAClone: "+ self.CDNAClone + " Ensemble.org: " + self.EnsemblOrg + \
		#	" LastChanceID: " + self.LastChanceID + \
		#	" " + self.DauerNumber1 + " " + self.DauerNumber2 + " " + self.DauerNumber3 + \
		#	" " + self.DauerNumber4 + "\n"

def updateIDStats(gene, stats):
	if gene.CDNAClone == "" and gene.WBGene == "" and gene.WB_Gene_ID == "" and gene.EnsemblOrg == "":
		stats['lastChanceId'] += 1
		stats['lastChanceTracking'].append(gene.LastChanceID)
	else:
		stats['reliable'] += 1
	stats['total'] += 1

def keyPrefixSearch(keysubstring, dict):
	for k in dict:
		if k.startswith(keysubstring):
			return k
	return None

##########################################################################
# First we read the sequence file and build up our identifier mapping

infile = open("C:\Users\Becky\Desktop\worm\c_elegans.current_development.cds_transcripts.fa","r")

line=infile.readline()
IDLookup = {}
while line != "":
	if line[0] == ">":
		line = line.strip()
		line = line.replace("gene=WBGene", "")
		line = line.replace(">", "")
		parts = line.split(" ")
		IDLookup[parts[0].lower()] = parts[1]
	line=infile.readline()
	
infile = open("C:\Users\Becky\Desktop\worm\c_elegans.current.cds_transcripts.fa","r")

line=infile.readline()
while line != "":
	if line[0] == ">":
		line = line.strip()
		line = line.replace("gene=WBGene", "")
		line = line.replace(">", "")
		parts = line.split(" ")
		IDLookup[parts[0].lower()] = parts[1]
	line=infile.readline()
	
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
			IDLookup[parts[1].lower()] = parts[0]
		if "not known" not in parts[2]:
			IDLookup[parts[2].lower()] = parts[0]
	line=infile.readline()

	
##########################################################################
# Now we process the data from dauerdb
			
currentStats = {}
currentStats['lastChanceId'] = 0
currentStats['lastChanceTracking'] = []
currentStats['reliable'] = 0
currentStats['total'] = 0

outfile=open("C:\Users\Becky\Desktop\worm\DauerOut", "w")

for file in os.listdir("C:\Users\Becky\Desktop\worm\dauer intruder"):
	if not isfile(join("C:\Users\Becky\Desktop\worm\dauer intruder", file)):
		continue
	if ".html" in file:
		continue
	if ".py" in file:
		continue
	
	infile = open(join("C:\Users\Becky\Desktop\worm\dauer intruder", file),"r")

	thisgene = Gene()
	ThisState = States.LookingFor_ViewDetailInfo

	line=infile.readline()
	while line != "":
		
		if ThisState == States.LookingFor_ViewDetailInfo:
			if "view_detail_info" in line:
				ThisState = States.LookingFor_ImgLongTermDauer1
			else:
				pass
				
		elif ThisState == States.LookingFor_ImgLongTermDauer1:
			if "img_long_term_dauer" in line:
				ThisState = States.LookingFor_ImgLongTermDauer2
			else:
				pass
				
		elif ThisState == States.LookingFor_ImgLongTermDauer2:
			if "img_long_term_dauer" in line:
				ThisState = States.LookingFor_ImgLongTermDauer3
			else:
				pass
				
		elif ThisState == States.LookingFor_ImgLongTermDauer3:
			if "img_long_term_dauer" in line:
				ThisState = States.LookingFor_ImgLongTermDauer4
			else:
				pass
				
		elif ThisState == States.LookingFor_ImgLongTermDauer4:
			if "img_long_term_dauer" in line:
				ThisState = States.LookingFor_SevenTabsTR1
			else:
				pass
				
		elif ThisState == States.LookingFor_SevenTabsTR1:
			if line.startswith("\t\t\t\t\t\t\t<tr>"):
				ThisState = States.LookingFor_SevenTabsTR2
			else:
				pass
				
		elif ThisState == States.LookingFor_SevenTabsTR2:
			if line.startswith("\t\t\t\t\t\t\t<tr>"):
				ThisState = States.LookingFor_SevenTabsTR3
			else:
				pass
				
		elif ThisState == States.LookingFor_SevenTabsTR3:
			if line.startswith("\t\t\t\t\t\t\t<tr>"):
				ThisState = States.LookingFor_SevenTabsTR4
			else:
				pass
				
		elif ThisState == States.LookingFor_SevenTabsTR4:
			if line.startswith("\t\t\t\t\t\t\t<tr>"):
				ThisState = States.LookingFor_DauerNumber1
			else:
				pass
				
		elif ThisState == States.LookingFor_DauerNumber1:
			if "<!--<br>(" not in line:
				raise Exception("I should have found a comment in this line")
			else:
				matches = re.search("\((-?\d\.?\d*)\)", line)
				if matches:
					thisgene.DauerNumber1 = matches.group(1)
					ThisState = States.LookingFor_DauerNumber2
				else:
					raise Exception("Should have gotten a regexp match")
				
		elif ThisState == States.LookingFor_DauerNumber2:
			if "<!--<br>(" in line:
				matches = re.search("\((-?\d\.?\d*)\)", line)
				if matches:
					thisgene.DauerNumber2 = matches.group(1)
					ThisState = States.LookingFor_DauerNumber3
				else:
					raise Exception("Should have gotten a regexp match")
					
		elif ThisState == States.LookingFor_DauerNumber3:
			if "<!--<br>(" in line:
				matches = re.search("\((-?\d\.?\d*)\)", line)
				if matches:
					thisgene.DauerNumber3 = matches.group(1)
					ThisState = States.LookingFor_DauerNumber4
				else:
					raise Exception("Should have gotten a regexp match")

		elif ThisState == States.LookingFor_DauerNumber4:
			if "<!--<br>(" in line:
				matches = re.search("\((-?\d\.?\d*)\)", line)
				if matches:
					thisgene.DauerNumber4 = matches.group(1)
					ThisState = States.LookingFor_Identifier
				else:
					raise Exception("Should have gotten a regexp match")
			
		elif ThisState == States.LookingFor_Identifier:
			if "WB_GENE_ID=" in line:
				matches = re.search("WB_GENE_ID=([^ ]+)", line)
				if matches:
					thisgene.WB_Gene_ID = matches.group(1).lower()
			if "cDNA clone" in line:
				matches = re.search("cDNA clone:?\s?([^ ]+)", line)
				if matches:
					thisgene.CDNAClone = matches.group(1).lower()
			if "WBGene" in line:
				matches = re.search("WBGene(\d+)", line)
				if matches:
					thisgene.WBGene = matches.group(1).lower()
			if "ensembl.org" in line:
				matches = re.search("\[<a href='http://www.ensembl.org/Caenorhabditis_elegans/geneview\?gene=([^']+)", line)
				if matches:
					thisgene.EnsemblOrg = matches.group(1).lower()
			if "REP_DB" in line:
				matches = re.search(">([^_ /]+)", line)
				if matches:
					thisgene.LastChanceID = matches.group(1).lower()
			
			if "view_detail_info" in line:
				if thisgene.hasIdentifier():
					if thisgene.MeetsConservativeReqs():
						outfile.write(thisgene.prettyOutput())
					updateIDStats(thisgene, currentStats)
					thisgene = Gene()
					ThisState = States.LookingFor_ImgLongTermDauer1
				else:
					raise Exception("What gene is this?")
			else:
				pass
				
		else:# ThisState gets corrupted somehow
			raise Exception("How did I get into this state?")
			
		line=infile.readline()

	if thisgene.MeetsConservativeReqs():
		outfile.write(thisgene.prettyOutput())
	updateIDStats(thisgene, currentStats)
	

print "Final Statistics"
print "Total Genes: ", currentStats['total']
print "Reliable Genes: ", currentStats['reliable']
print "Last Chance IDs: ", currentStats['lastChanceId']
for id in currentStats['lastChanceTracking']:
	print "\t", id