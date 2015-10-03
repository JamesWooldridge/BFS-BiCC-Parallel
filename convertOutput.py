files = [
	"100-10"#, "100-50",
	# "1000-10", "1000-250",
	# "10000-50", "10000-250",
	# "100000-200", "100000-25000",
	# "1000000-250", "1000000-250000"
]

for filename in files:
	inputFilename = "output/" + filename + "-par.bic"
	outputFilename = "output/c_" + filename + "-par.bic"

	components = {}
	with open(inputFilename, "r") as inputFile:
		for line in inputFile:
			line = line.rstrip()
			line = line.split(" = ")
			label = line[1]
			if not label in components:
				components[label] = []
			
			components[label].append(line[0])

	with open(outputFilename, "w") as outputFile:
		for value in components.itervalues():
			for edge in value:
				outputFile.write(edge + "\n")
			outputFile.write("\n")