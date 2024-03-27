#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

class GeneSequencing:

	def __init__( self ):
		pass

	
# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length
		results = []

		for i in range(len(sequences)):
			jresults = []
			for j in range(len(sequences)):

				if(j < i):
					s = {}
				else:
###################################################################################################
# your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
					if(not banded): score = self.unbandedAlignment(align_length, sequences[i], sequences[j])
					
					else: score = self.bandedAlignment(align_length, sequences[i], sequences[j])
					
					alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
						len(sequences[i]), align_length, ',BANDED' if banded else '')
					alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
						len(sequences[j]), align_length, ',BANDED' if banded else '')
###################################################################################################					
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
					table.update()	
				jresults.append(s)
			results.append(jresults)
		return results
	

	def unbandedAlignment( self, align_length, sequence1, sequence2):
		if( len(sequence1) < align_length  and  len(sequence2) < align_length):	# if both sequence 1 and sequence 2 are smaller than the alignment length
			table = []
			for i in range(len(sequence1)+1):
				newArray = []
				for j in range(len(sequence2)+1):
					newArray.append(0)
				table.append(newArray)

			for i in range(len(sequence1)+1):
				table[i][0] = i
			for j in range(len(sequence2)+1):
				table[0][j] = j

			for i in range(1, len(sequence1)+1):
				for j in range(1, len(sequence2)+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
					'''
					if(score1 <= score2 and score1 <= score3):
						table[i][j] = score1
					elif(score2 < score1 and score2 <= score3):
						table[i][j] = score2
					else: table[i][j] = score3
					'''
			return table[len(sequence1)][len(sequence2)]
		
		elif(len(sequence1) < align_length and len(sequence2) >= align_length):	# if sequence 1 is less than the align length
			table = []
			for i in range(len(sequence1)+1):
				newArray = []
				for j in range(align_length+1):
					newArray.append(0)
				table.append(newArray)

			for i in range(len(sequence1) + 1):
				table[i][0] = i
			for j in range(align_length + 1):
				table[0][j] = j

			for i in range(1, len(sequence1)+1):
				for j in range(1, align_length+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
					'''
					if(score1 <= score2 and score1 <= score3):
						table[i][j] = score1
					elif(score2 < score1 and score2 <= score3):
						table[i][j] = score2
					else: table[i][j] = score3
					'''
			return table[len(sequence1)][align_length]
		
		elif( len(sequence1) >= align_length  and  len(sequence2) < align_length):	# if sequence 2 is less than the align length
			table = []
			for i in range(align_length + 1):
				newArray = []
				for j in range(len(sequence2)+1):
					newArray.append(0)
				table.append(newArray)

			for i in range(align_length+1):
				table[i][0] = i
			for j in range(len(sequence2)+1):
				table[0][j] = j

			for i in range(1, align_length+1):
				for j in range(1, len(sequence2)+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
					'''
					if(score1 <= score2 and score1 <= score3):
						table[i][j] = score1
					elif(score2 < score1 and score2 <= score3):
						table[i][j] = score2
					else: table[i][j] = score3
					'''
			return table[align_length][len(sequence2)]


		else:	# if both are bigger than the align length
			table = []
			for i in range(align_length+1):
				newArray = []
				for j in range(align_length+1):
					newArray.append(0)
				table.append(newArray)

			for i in range(align_length+1):
				table[i][0] = i
			for j in range(align_length+1):
				table[0][j] = j

			for i in range(1, align_length+1):
				for j in range(1, align_length+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					'''
					if(score1 <= score2 and score1 <= score3):
						table[i][j] = score1
					elif(score2 < score1 and score2 <= score3):
						table[i][j] = score2
					else: table[i][j] = score3
					'''
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
			return table[align_length][align_length]
		
	
	def bandedAlignment( self, align_length, sequence1, sequence2):
		return



