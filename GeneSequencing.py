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
			# This section initializes the chart being filled out
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

			# Compares the 2 sequences and assigns them different scores
			for i in range(1, len(sequence1)+1):
				for j in range(1, len(sequence2)+1):
					'''
					MATCH/SUB = score1/diagonal
					INSERT = score3/Left
					DELETE = score2/Top
					'''
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					# Determines which score is lowest and handles the tie breakers
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
					
			return table[len(sequence1)][len(sequence2)]
		
		elif(len(sequence1) < align_length and len(sequence2) >= align_length):	# if sequence 1 is less than the align length
			# This section initializes the chart being filled out
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

			# Compares the 2 sequences and assigns them different scores
			for i in range(1, len(sequence1)+1):
				for j in range(1, align_length+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					# Determines which score is lowest and handles the tie breakers
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
			return table[len(sequence1)][align_length]
		
		elif( len(sequence1) >= align_length  and  len(sequence2) < align_length):	# if sequence 2 is less than the align length
			# This section initializes the chart being filled out
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

			# Compares the 2 sequences and assigns them different scores
			for i in range(1, align_length+1):
				for j in range(1, len(sequence2)+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					# Determines which score is lowest and handles the tie breakers
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
			return table[align_length][len(sequence2)]


		else:	# if both are bigger than the align length
			# This section initializes the chart being filled out
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

			# Compares the 2 sequences and assigns them different scores
			for i in range(1, align_length+1):
				for j in range(1, align_length+1):
					if(sequence1[i-1] == sequence2[j-1]):
						score1 = table[i-1][j-1] + MATCH
					else: score1 = table[i-1][j-1] + SUB
					score2 = table[i][j-1] + INDEL
					score3 = table[i-1][j] + INDEL
					# Determines which score is lowest and handles the tie breakers
					if(score3 <= score1 and score3 <= score2):
						table[i][j] = score3
					elif(score2 <= score1 and score2 < score3):
						table[i][j] = score2
					else: table[i][j] = score1
			return table[align_length][align_length]
		
	
	def bandedAlignment( self, align_length, sequence1, sequence2):
		arrayWidth = 2*MAXINDELS + 1
		table = []
		if(len(sequence1)<=align_length):
			for first_word_index in range(len(sequence1)):
				newArray = []
				for second_word_index in range(arrayWidth):
					if(first_word_index < MAXINDELS  and  second_word_index < MAXINDELS - first_word_index):
						newArray.append(float('inf'))
					elif(first_word_index > 2*MAXINDELS  and  second_word_index >= arrayWidth - MAXINDELS and second_word_index > len(sequence1) - first_word_index + MAXINDELS -1): 
						newArray.append(float('inf'))
					else: newArray.append(0)
				table.append(newArray)

			for second_word_index in range(arrayWidth):
				if(table[0][second_word_index] == float('inf')): continue
				else: table[0][second_word_index] = second_word_index - MAXINDELS

			second_word_index = MAXINDELS
			for first_word_index in range(MAXINDELS + 1):
				table[first_word_index][second_word_index] = first_word_index
				second_word_index -= 1

			for first_word_index in range(1, len(sequence1)):
				if(first_word_index - MAXINDELS <= 0): second_word_index_string = 0
				else: second_word_index_string = first_word_index - MAXINDELS
				for second_word_index in range(arrayWidth):
					if(table[first_word_index][second_word_index] == float('inf')): continue	# This means the value is not included in our band
					elif(table[first_word_index][second_word_index - 1] == float('inf')): continue	# This means the value is on the upper left part of the band, so we don't update its value
					else:
						if((second_word_index == 0 and table[first_word_index][second_word_index] != float('inf')) or (second_word_index == MAXINDELS - first_word_index)):	# This means that the value is on the left side of the band
							if(sequence1[first_word_index] == sequence2[second_word_index_string]):
								match_sub_score = table[first_word_index-1][second_word_index] + MATCH
							else:
								match_sub_score = table[first_word_index-1][second_word_index] + SUB
							delete_score = table[first_word_index-1][second_word_index+1] + INDEL
							if(delete_score <= match_sub_score): table[first_word_index][second_word_index] = delete_score
							else: table[first_word_index][second_word_index] = match_sub_score
						elif((second_word_index == arrayWidth-1) or (table[first_word_index][second_word_index + 1] == float('inf'))):	#This means that the value is on the right side of the band.
							if(sequence1[first_word_index] == sequence2[second_word_index_string]):
								match_sub_score = table[first_word_index-1][second_word_index] + MATCH
							else:
								match_sub_score = table[first_word_index-1][second_word_index] + SUB
							insert_score = table[first_word_index][second_word_index -1] + INDEL
							if(insert_score <= match_sub_score): table[first_word_index][second_word_index] = insert_score
							else: table[first_word_index][second_word_index] = match_sub_score
						else:
							if(sequence1[first_word_index] == sequence2[second_word_index_string]):
								match_sub_score = table[first_word_index-1][second_word_index] + MATCH
							else:
								match_sub_score = table[first_word_index-1][second_word_index] + SUB
							delete_score = table[first_word_index-1][second_word_index+1] + INDEL
							insert_score = table[first_word_index][second_word_index -1] + INDEL
							if(insert_score <= match_sub_score and insert_score <= delete_score):
								table[first_word_index][second_word_index] = insert_score
							if(delete_score <= match_sub_score):
								table[first_word_index][second_word_index] = delete_score
							else:
								table[first_word_index][second_word_index] = match_sub_score
						second_word_index_string += 1
			print(table)
			print('\r\n')
			return table[len(sequence1)-1][MAXINDELS]
		
		else:
			#table = []
			for first_word_index in range(align_length):
				newArray = []
				for second_word_index in range(arrayWidth):
					if(first_word_index < MAXINDELS  and  second_word_index < MAXINDELS - first_word_index):
						newArray.append(float('inf'))
					elif(first_word_index > 2*MAXINDELS  and  second_word_index >= arrayWidth - MAXINDELS and second_word_index > align_length - first_word_index + MAXINDELS -1): 
						newArray.append(float('inf'))
					else: newArray.append(0)
				table.append(newArray)
			
			for second_word_index in range(arrayWidth):
				if(table[0][second_word_index] == float('inf')): continue
				else: table[0][second_word_index] = second_word_index - MAXINDELS

			second_word_index = MAXINDELS
			for first_word_index in range(MAXINDELS + 1):
				table[first_word_index][second_word_index] = first_word_index
				second_word_index -= 1
			
			for first_word_index in range(1, align_length):
				if(first_word_index - MAXINDELS <= 0): second_word_index_string = 0
				else: second_word_index_string = first_word_index - MAXINDELS
				for second_word_index in range(arrayWidth):
					if(table[first_word_index][second_word_index] == float('inf')): continue	# This means the value is not included in our band
					elif(table[first_word_index][second_word_index - 1] == float('inf')): continue	# This means the value is on the upper left part of the band, so we don't update its value
					else:
						if((second_word_index == 0 and table[first_word_index][second_word_index] != float('inf')) or (second_word_index == MAXINDELS - first_word_index)):	# This means that the value is on the left side of the band
							if(sequence1[first_word_index] == sequence2[second_word_index_string]):
								match_sub_score = table[first_word_index-1][second_word_index] + MATCH
							else:
								match_sub_score = table[first_word_index-1][second_word_index] + SUB
							delete_score = table[first_word_index-1][second_word_index+1] + INDEL
							if(delete_score <= match_sub_score): table[first_word_index][second_word_index] = delete_score
							else: table[first_word_index][second_word_index] = match_sub_score
						elif((second_word_index == arrayWidth-1) or (table[first_word_index][second_word_index + 1] == float('inf'))):	#This means that the value is on the right side of the band.
							if(sequence1[first_word_index] == sequence2[second_word_index_string]):
								match_sub_score = table[first_word_index-1][second_word_index] + MATCH
							else:
								match_sub_score = table[first_word_index-1][second_word_index] + SUB
							insert_score = table[first_word_index][second_word_index -1] + INDEL
							if(insert_score <= match_sub_score): table[first_word_index][second_word_index] = insert_score
							else: table[first_word_index][second_word_index] = match_sub_score
						else:
							if(sequence1[first_word_index] == sequence2[second_word_index_string]):
								match_sub_score = table[first_word_index-1][second_word_index] + MATCH
							else:
								match_sub_score = table[first_word_index-1][second_word_index] + SUB
							delete_score = table[first_word_index-1][second_word_index+1] + INDEL
							insert_score = table[first_word_index][second_word_index -1] + INDEL
							if(insert_score <= match_sub_score and insert_score <= delete_score):
								table[first_word_index][second_word_index] = insert_score
							if(delete_score <= match_sub_score):
								table[first_word_index][second_word_index] = delete_score
							else:
								table[first_word_index][second_word_index] = match_sub_score
					second_word_index_string += 1
			return table[align_length-1][arrayWidth-MAXINDELS-1]



