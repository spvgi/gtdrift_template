import sys

# This script is used to generate the json and query used in json files of other scripts
# python3 generate_json_and_query.py data/resources/organisms_list conf_assembly_list conf_query 

dico = {}
with open(sys.argv[1]) as reader:
    data = reader.readlines()
    for line in data[1:]:
        acc = line.split(".")[0][4:]
        if acc in dico:
            dico[acc].append(line.rstrip())
        else :
            dico[acc] = [line.rstrip()]
unic_dico = {}
for key in dico:
    val  = dico[key]
    if len(val) == 1 :
        unic_dico[key] = val[0]
    if len(val) == 2 :
        test = 0
        for candidat in val:
            if candidat[0:3] == "GCF" :
                unic_dico[key] = candidat
                test += 1
        if test != 1 :
            sys.exit("we want one and only one  GCF")
    if len(val) > 2 :
        sys.exit("too many assemblies")
accessions = list(unic_dico.values())
with open(sys.argv[2], 'w') as writer, open(sys.argv[2]+".col", 'w') as writer_col, open(sys.argv[3], 'w') as writer2, open(sys.argv[3]+".col", 'w') as writer2_col:
    writer.write("[ ")
    writer.write('"'+accessions[0]+'" ')
    writer_col.write("[ ")
    writer_col.write('"'+accessions[0]+'" ')
    writer2.write('('+accessions[0].split(".")[0]+') ')
    writer2_col.write('('+accessions[0].split(".")[0]+') ')
    for line in accessions[1:]:
            writer.write(', "'+line+'" ')
            writer_col.write(', \n"'+line+'" ')
            writer2.write('OR ('+line.split(".")[0]+') ')
            writer2_col.write('\n OR ('+line.split(".")[0]+') ')
    writer.write("]")
    writer_col.write("]")
