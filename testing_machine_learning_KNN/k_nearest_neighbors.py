import sys
import csv
import random
import math
import operator


def euclideanDistance(instace1, instance2, length):
	distance = 0
	for x in range(length):
		distance += pow(float(instace1[x]) - float(instance2[x]), 2)
	return math.sqrt(distance)

def get_neighbors(training_set, test_set, k):
	distances = []
	length = len(test_set) - 1
	for x in range(len(training_set)):
		dist = euclideanDistance(test_set, training_set[x], length)
		distances.append((training_set[x], dist))
	#distances.sort(key=operator.itemgetter(1), reverse=True) -> En caso de comprobar los lejanos 
	distances.sort(key=operator.itemgetter(1), reverse=False)
	neighbors = []
	for x in range(k):
		neighbors.append(distances[x][0][-1])
	return neighbors

def getResponse(neighbors):
	classVotes = {}
	for x in range(len(neighbors)):
		response = neighbors[x][-1]
		if response in classVotes:
			classVotes[response] += 1
		else:
			classVotes[response] = 1
	sortedVotes = sorted(classVotes.iteritems(), key=operator.itemgetter(1), reverse=True)
	return sortedVotes[0][0]
 
def getAccuracy(testSet, predictions):
	correct = 0
	for x in range(len(testSet)):
		if testSet[x][-1] == predictions[x]:
			correct += 1
	return (correct/float(len(testSet))) * 100.0

def loadDataset(lines, split, trainingSet=[] , testSet=[]):

    dataset = list(lines)
    for x in range(len(dataset)-1):
        for y in range(4):
            dataset[x][y] = float(dataset[x][y])
        if random.random() < split:
            trainingSet.append(dataset[x])
        else:
            testSet.append(dataset[x])


label_dict = {}
with open('test_labels.txt') as f:
	lines = f.readlines()
	for line in lines:
		name, label = line.split(',')
		label_dict[name] = label

training_set = []
with open('outmash_training_59_labels.csv', 'rb') as csvfile:
	lines = csv.reader(csvfile)
	for row in lines:
		values, name = row[:-1], row[-1]
		try:
			label = label_dict[name].strip()
		except KeyError:
			label = 'i'

		training = values
		training.insert(len(training), label)
		training_set.append(training)

def main():

	trainingSet=[]
	testSet=[]
	split = 0.67
	loadDataset(training_set, split, trainingSet, testSet)
	print
	print('Train set: ' + repr(len(trainingSet)))
	print('Test set: ' + repr(len(testSet)))
	# generate predictions
	predictions=[]
	k = 9
	print
	for x in range(len(testSet)):
		neighbors = get_neighbors(trainingSet, testSet[x], k)
		result = getResponse(neighbors)
		predictions.append(result)
		print('> predicted=' + repr(result) + ', actual=' + repr(testSet[x][-1]))
	accuracy = getAccuracy(testSet, predictions)
	print('\n')
	print('k-fold cross-validation: ', k)
	print('Accuracy: ' + repr(accuracy) + '%')

if __name__ == '__main__':
	main()


