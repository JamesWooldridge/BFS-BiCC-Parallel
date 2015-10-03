from sets import Set

files = [
	"100-10", "100-50",
	"1000-10", "1000-250",
	"10000-50", "10000-250",
	"100000-200", "100000-25000",
	"1000000-250", "1000000-250000"
]

for filename in files:
	correctFilename = "output/" + filename + ".bic"
	testFilename = "output/" +  filename + "-par.bic"

	correctComponents = Set()
	currSet = Set()
	with open(correctFilename, "r") as correctFile:
		for line in correctFile:
			if line == "\n":
				correctComponents.add(currSet)
				currSet = Set()
			else:
				currSet.add(line.rstrip())

	foundComponents = Set()
	currSet = Set()
	with open(testFilename, "r") as testFile:
		for line in testFile:
			if line == "\n":
				foundComponents.add(currSet)
				currSet = Set()
			else:
				currSet.add(line.rstrip())

	print correctComponents == foundComponents