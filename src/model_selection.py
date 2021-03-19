import argparse
import pandas as pd
import numpy as np

def get_arguments():
    """
    Lecture des arguments
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-p', '--path', help="Path to HHR file")
    parser.add_argument('-f', '--fasta', help = "Path to the fasta file analysed")
    args = parser.parse_args()
    return args

def get_stats(raw_stat, name, length):
	"""
	Fonction extrayant les statistiques proposé pour un alignement HHR et les renvoyants sous forme de liste
	----------
	Input :
		raw_stat = ligne de statistique brut en chaine de caractère
		name = Nom de la séquence
		length = Longeur de la séquence d'entrée
	----------
	Output :
		stats = liste de valeurs des statistiques obtenus
	----------
	"""
	stats=[]
	stats.append(name) #Nom
	stats.append(float(raw_stat[1])) #Prob
	stats.append(float(raw_stat[3])) #E-values
	stats.append(float(raw_stat[5])) #Score
	stats.append(int(raw_stat[7])/length*100) #Couvrance
	stats.append(float(raw_stat[9].strip("%"))) #Identité
	stats.append(float(raw_stat[11])) #Similarité
	return stats


def pdb_align_parser(raw_align, length):
	"""
	Function use to transform a raw_pdb_alignment in an understable format to evaluate
	----------
	Input :
		raw_stat = ligne de statistique brut en chaine de caractère
		length = Longeur de la séquence d'entrée
	----------
	Output :
		stats = liste de valeurs des statistiques obtenus, contennat les prédictions de structures secondaires
	----------
	"""
	name = raw_align[0][1:7]
	stats = get_stats(raw_align[1].replace('=', ' ').split(), name, length)
	#Stating position of the alignment
	debut = int(raw_align[3].split()[2])
	stats.append(debut)
	seq = []
	dssp = []
	ss_pred = []
	ss_seqali = []
	for i, line in enumerate(raw_align):
		if line.startswith("Q sp|"):
			seq += line.split()[3]
		flag = 5
		if line.startswith("T ss_dssp"):
			dssp += line.split()[2]
			flag = 4
		if line.startswith("T ss_pred"):
			ss_pred += line.split()[2]
			ss_seqali += raw_align[i-flag].split()[3]

	for i in range(len(seq)):
		if seq[i] == "-":
			ss_pred[i] = "."
	ss_pred = ''.join(ss_pred)
	fin_dssp = debut + len(dssp)
	stats.append(fin_dssp)
	#Cleaning insertion
	ss_pred = ''.join(ss_pred.split("."))
	fin_pred = debut + len(ss_pred)
	stats.append(fin_pred)
	stats.append(dssp)
	stats.append(ss_pred)
	stats.append(ss_seqali)
	return stats

def itasser_ss_file(ssList, fasta_seq, file_name):
	"""
	Fonction écrivant dans file_name les structures secondaires attendus pour une séquence fasta
	----------
	Input :
		ssList = liste des structures secondaires
		fasta_seq = Séquence fasta
		file_name = fichier de sortie
	----------
	"""
	with open(file_name, "w") as fillout:
		for i, (aa, ss) in enumerate(zip(fasta_seq, ssList)):
			if ss != " ":
				if ss == "e" or ss == "E":
					fillout.write("{:<6d} {:<} {:<}\n".format(i+1, aa, "S"))
				else:
					fillout.write("{:<6d} {:<} {:<}\n".format(i+1, aa, ss.upper()))


def main():
	args = get_arguments()
	path = args.path
	fasta = args.fasta
	fasta_name = fasta.split(".")[0]
	print(fasta, path)

	with open(fasta, "r") as f:
		fasta_seq = f.read()
	sequence = "".join(fasta_seq.split("\n")[1:])

	#Save the hhr file in RAM memmory
	with open(path, "r") as f:
		full_file = f.readlines()

	#Get sequence length
	length = int(full_file[1].split()[1])
	print("longeur de {}".format(length))

	# Transformer le hhr en un format lisible
	df = pd.DataFrame(columns=["Name", "Prob", "Eval", "Score", "Couvrance", "Identities", "Similarity", "debut", "fin_dssp", "fin_pred", "dssp", "ss_pred", "ss_seqalign"])
	flag=0
	aliList=[]
	alignement=[]
	for line in full_file:
		#New_align
		if line.startswith(">"):
			#End of one align, add them to global lsit and reset
			if flag==1:
				aliList.append(alignement)
				alignement=[]
				alignement.append(line)
			#Skip header
			if flag==0:
				flag=1
				alignement.append(line)
		#Main body of align
		elif flag == 1:
			alignement.append(line)

	for i, line in enumerate(aliList):
		df.loc[i]=pdb_align_parser(line, length)

	""" ########### #######   ######## ########    ##  ###    ##   	########
			##      ##    ##  ##       ##     ##   ##  ## #  ##    ##
			##      #######   ####     ##      ##  ##  ##  # ##    ##  #######
			##      ##   ##   ##       ##     ##   ##  ##   ###    ##    ##
			##      ##    ##  ######## ########    ##  ##    ##     ######  
	"""
	ss_model = [" "] * length #Create an empty liste with define length
	name = []
	good_flag = 0
	for index, row in df.iterrows():
		#Check si la strcuture est suffisament bonne pour être utilisée seule en guise de modèle
		if row['Couvrance'] > 80 and row['Identities'] > 40 and row['Eval'] < 1e-10:
			name = row['Name']
			for i, x in enumerate(range(row['debut'],row['fin_pred'])):
				if row['ss_pred'][i] != "-":
					ss_model[x]=row['ss_pred'][i]
			itasser_ss_file(ss_model, sequence, row['Name']+"_ss.txt")
			good_flag = 1
			break
		#Multiple
		#Sinon on va construire le modèle en forceant les structures secondaires a partir d'alignements multiples
		if row['Couvrance'] > 20 and row['Identities'] > 10 and row['Eval'] < 1e-10:
			good_flag = 2
			couv = 0
			longest_empty = 0
			start_empty = 0
			ss_pred = []
			#Calcul de la plus long portion de l'alignement pas déja couverte par un alignement de meilleurqualitée
			for i in range(row['debut'], row['fin_pred']):
				#Si l'alignement s'intéresse a une partie pas encore comblé, on enregiste
				if ss_model[i] != "":
					couv += 1
				#On garde la plus grande aprtie vide de l'alignement
				elif couv != 0 and couv > longest_empty:
					longest_empty = couv
					start_empty = i-couv
					couv = 0
			#Si la partie vite est l'ensemble de l'alignement alors ajout
			if couv != 0 and longest_empty == 0:
				name.append(row['Name'])
				longest_empty = couv
				start_empty = row['debut']

			#Ajout a la séquence de structure secondaire prédite l'alignement calculé
			elif longest_empty > 20 :  #Critère de qualité
				ss_pred = row['ss_pred'][start_empty-row['debut']:start_empty+longest_empty-row['debut']]
				name.append(row['Name'])
				#Ajout dans le modèle le fragment trouvé
				for i,j in enumerate(range(start_empty, start_empty+longest_empty-2)):
					print(i, len(ss_pred))
					if ss_pred[i] != "-":
						ss_model[j] = ss_pred[i].upper()


	if good_flag == 1:
		print("cas optimal lancer i-tasser avec le fasta {} et la séquence modèle {}".format(fasta, name.split("_")[0]))
	if good_flag == 2:
		fillout = path.split(".")[0] + "_SS.txt"
		print("Cas non optimal,lancer itasser avec fasta : {}\t et modèle de structure secondaire: {}".format(fasta, fillout))
		print("Les alignements utilisé pour al construction du modèle sont : {}".format(name))
		itasser_ss_file(ss_model, sequence, fillout)
	else :
		return 0

if __name__ == "__main__":
	main()

