import sys

# This script is used to generate the json and query used in json files of other scripts

with open(sys.argv[1]) as reader, open(sys.argv[2], 'w') as writer, open(sys.argv[3], 'w') as writer2:
    data = reader.readlines()
    writer.write("[ ")
    writer.write('"'+data[1].rstrip()+'" ')
    writer2.write('('+data[1].split(".")[0]+') ')    
    for line in data[2:]:
            writer.write(', "'+line.rstrip()+'" ')
            writer2.write('OR ('+line.split(".")[0]+') ')
    writer.write("]")
