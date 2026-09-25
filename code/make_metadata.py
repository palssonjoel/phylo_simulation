
files = ["../output/simulation/1/subsamples/random_150_no1.fasta", "../output/simulation/1/subsamples/random_500_no1.fasta"]
headers = []

for file in files:
  outfile = (file.split("/")[5]).split(".")[0] + "_metadata.csv"

  with open(file, "r") as file, open(outfile, "w") as outfile:
    outfile.write("name,location,time\n")
    for line in file:
      if line.startswith(">"):
        line = line.replace(">", "")

        split_line = line.split("|")
        limit = len(split_line)
        n = 0

        for i in split_line:
          
          if n < limit-1:
            i = i.strip()
            outfile.write(f"{i},")
          else: 
            i = i.strip()
            outfile.write(f"{i}\n")
          n += 1


