#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 07:50:14 2018
adopted from orf-finder.py
"""
import sys 
import re
import os
from Bio import SeqIO


fileName = sys.argv[1]
outputFileName = sys.argv[2]
minProLen = 32
maxProLen = 100
hitsToPull =3
altCodonStringency = 49
noCodonStringency = 99
sequenceType = 'both'
replace = 'n'
force = 'n'
unresolvedCodon = 0
THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
my_file = os.path.join(THIS_FOLDER, fileName)



xRegex = re.compile(r'X+')   


if outputFileName != None:
        if sequenceType.lower() != 'both':
                if os.path.isfile(outputFileName) and force.lower() != 'y':
                        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete this older file, or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(outputFileName) and force.lower() == 'y':
                        os.remove(outputFileName)
        else:
                outPrefix = outputFileName.rsplit('.',1)
                protOutName = outPrefix[0] + '_prot.' + outPrefix[1]
                nuclOutName = outPrefix[0] + '_nucl.' + outPrefix[1]
                if os.path.isfile(protOutName) and force.lower() != 'y':
                        print('There is already a file named ' + protOutName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(protOutName) and force.lower() == 'y':
                        os.remove(protOutName)
                if os.path.isfile(nuclOutName) and force.lower() != 'y':
                        print('There is already a file named ' + nuclOutName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(nuclOutName) and force.lower() == 'y':
                        os.remove(nuclOutName)

if unresolvedCodon != 0:
        print('Program has noted that you are allowing the discovery of ORFs with unresolved codon regions. This is risky behaviour, since this program cannot guarantee that an unresolved region does not contain a stop codon. Subsequently, you can have chimeras form from two separate ORFs. YOU MUST VERIFY ANY ORFS WITH UNRESOLVED REGIONS! The best way to do this is with BLAST against homologous proteins. You have been warned.')

# Check for silly settings
if hitsToPull == 0:
        print('You set numhits to 0. There\'s no point running this program if you aren\'t getting any output! Specify a number >0 and try again.')
        quit()

# Load the fasta file as a generator object, get the total number of sequences in the file, then re-load it for the upcoming loop
input_handle = open(my_file,'rU')
records = SeqIO.parse(input_handle, "fasta")
totalCount = 0
for record in records:
        totalCount += 1
        
input_handle = open(my_file,'rU')    
records = SeqIO.parse(input_handle, 'fasta')

### CORE PROCESSING LOOP
print('Starting the core processing of this script now. Progress bar is displayed below. Every 10,000 sequences, current progress will be saved to the output file(s) to reduce memory footprint.')

# Declare overall values needed before loop start
startCodon = re.compile(r'^.*?(M.*)')           # Regex to pull out just the sequence starting with a Methionine (or ATG)
ongoingCount = 0
outputProt = []                                 # These get reset whenever we output to file
outputNucl = []
rememberPrint = -1

# Get the nucleotide (record) out of our generator (records) and grab them ORFs!
for record in records:
        # Progress bar
        progress = ((ongoingCount+1)/totalCount)*100
        if int(progress)%1==0 and rememberPrint != int(progress):
            print('|' + str(int(progress)) + '% progress|' + str(ongoingCount+1) + ' sequences scanned for ORFs'+'\r')      
            # Need to +1 to ongoingCount to counteract 0-index
            rememberPrint = int(progress) 
        # Declare output holding values that should reset for each transcript/record
        tempOverallProt = []
        tempOverallNucl = []
        tempMProt = []
        tempMNucl = []
        tempAltProt = []
        tempAltNucl = []
        tempNoneProt = []
        tempNoneNucl = []
        # Parental loop
        for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                for frame in range(3):
                        length = 3 * ((len(record)-frame) // 3)
                        frameNuc = str(nuc[frame:frame+length])
                        frameProt = str(nuc[frame:frame+length].translate(table=1))
                        # Split protein/nucleotide into corresponding ORFs
                        ongoingLength = 0                                       # The ongoingLength will track where we are along the unresolvedProt sequence for getting the nucleotide sequence
                        splitNucleotide = []
                        splitProtein = []
                        frameProt = frameProt.split('*')
                        for i in range(len(frameProt)):
                                if len(frameProt) == 1 or i + 1 == len(frameProt):    # This means the splitProtein has no stop codons or we're looking at the last ORF which also has no stop codon
                                        splitProtein.append(frameProt[i])
                                        splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i])*3])      # This will grab the corresponding nucleotide region
                                        ongoingLength += len(frameProt[i])*3
                                else:
                                        splitProtein.append(frameProt[i] + '*')
                                        splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i] + '*')*3])       
                                        ongoingLength += (len(frameProt[i]) + 1)*3
                        
                        # Fix unresolved regions
                        resolvedProt = []
                        resolvedNuc = []
                        indicesForDel = []                                                              # We'll hold onto the indices of any splitProtein components that have X's in them. We'll loop through this later to delete them from splitProtein/Nucleotide, then we'll add the resolved segments to splitProtein/Nucleotide
                        for i in range(len(splitProtein)):
                                if 'X' in splitProtein[i]:
                                        posProt = []
                                        for x in re.finditer(xRegex, splitProtein[i]):
                                                if x.end() - x.start() > unresolvedCodon:
                                                        posProt += [x.start(), x.end()]
                                        if posProt == []:                                               # If posProt still == [], that means we didn't find any unresolved regions that exceed our cut-off
                                                continue
                                        indicesForDel.insert(0, i)                                      # Insert it at 0 so we don't need to sort it at the end [we need to loop through a reversed list so we can directly delete the indices without messing up the order of splitProtein/Nucleotide]
                                        # Pull out resolved regions
                                        resolvedProt.append(splitProtein[i][:posProt[0]])               # We loop through our posProt by first grabbing everything before our first unresolved region
                                        resolvedNuc.append(splitNucleotide[i][:posProt[0]*3])
                                        for x in range(1, len(posProt)-1, 2):                           # We now, by skipping the first and last entry in posProt, can compare every coordinate pair that corresponds to a resolved region
                                                start = posProt[x]
                                                end = posProt[x+1]
                                                resolvedProt.append(splitProtein[i][start:end])
                                                resolvedNuc.append(splitNucleotide[i][start*3:end*3])
                                        resolvedProt.append(splitProtein[i][posProt[-1]:])              # We can now grab everything after our last unresolved region. If there was only one unresolved region, we never enter the above loop and just use the coordinate pair to get our start and end sequences
                                        resolvedNuc.append(splitNucleotide[i][posProt[-1]*3:])
                        # Delete old entries and add resolved entries
                        for index in indicesForDel:
                                del splitProtein[index]
                                del splitNucleotide[index]
                        splitProtein += resolvedProt                                                    # If we don't find any unresolved regions we wanted to delete, resolvedProt will be empty so nothing happens
                        splitNucleotide += resolvedNuc

                        # Enter the main processing loop with our resolved regions
                        for i in range(len(splitProtein)):                                                      # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                                # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                                mPro = ''
                                altPro = ''
                                nonePro = ''
                                topHit = ''
                                codonIndex = None
                                noneCodonContingency = None
                                # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                                if len(splitProtein[i]) < minProLen:                    # Disregard sequences that won't meet the size requirement without further processing
                                        continue
                                elif maxProLen != 0 and len(splitProtein[i]) > maxProLen:
                                        continue
                                acceptedPro = str(splitProtein[i])
                                # Alternative start coding      
                                nucSeqOfProt = splitNucleotide[i]                       # Don't need to do it, but old version of script extensively uses this value and cbf changing it
                                codons = re.findall('..?.?', nucSeqOfProt)              # Pulls out a list of codons from the nucleotide
                                for codon in codons:                                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                                        if codon == 'GTG' or codon == 'TTG':
                                                codonIndex = codons.index(codon)        # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                                break
                                        elif codon == 'CTG':
                                                if noneCodonContingency == None:        # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                                                        noneCodonContingency = codons.index(codon)

                                # Get the three ORF versions from each region inbetween stop codons
                                if 'M' in str(acceptedPro):                             # Obtains a traditional methionine initiated ORF starting from the first methionine if there is one in the sequence
                                        mPro = startCodon.search(str(acceptedPro)).groups()[0]  # Note that startCodon was declared at the start of this file         

                                if codonIndex != None:                                  # Gets the start position of the protein if we found a likely alternative start (aka a 'GTG' or 'TTG')
                                        altPro = acceptedPro[codonIndex:]
                                        if replace.lower() == 'y':
                                                altPro = 'M' + altPro[1:]               # If the argument is provided, this script assumes that the alternative start will be substituted with a methionine post-transcription
                                elif noneCodonContingency != None:                      # This will match an alternative start to 'CTG' only if 'TTG' or 'GTG' are not present
                                        altPro = acceptedPro[codonIndex:]
                                        if replace.lower() == 'y':
                                                altPro = 'M' + altPro[1:]

                                if i == 0:                                              # nonePro makes an assumption that the start of the ORF was not assembled properly resulting in the real start codon being cut off. Our stringency values will assess the likelihood of this hypothesis.
                                        nonePro = acceptedPro                           # Additionally, by only obtaining a 'nonePro' when it is in a protein fragment at the start of a frame (i.e., splitProtein[0]), we also make the (reasonable) assumption that any ORF inbetween two stop codons should itself have a start codon. This doesn't always hold true due to transcript assembly errors, but it must be assumed for the purpose of this script.                          

                                # Pull out the top hit from this protein fragment based upon how strict we want to be with accepting a traditional, alternative, or no codon start                                
                                if len(nonePro) > len(altPro) + noCodonStringency and len(nonePro) > len(mPro) + noCodonStringency:             # By adding on the stringency values declared earlier to the length of the protein, we can determine whether we want to consider an ORF without a start codon as legitimate
                                        topHit = nonePro
                                elif len(altPro) > len(mPro) + altCodonStringency:                                                              # Adding the stringency values here allows us to determine whether we can increase ORF length significantly by assuming an alternative start rather than a methionine start
                                        topHit = altPro
                                else:
                                        topHit = mPro                                                                                           # This is the default position unless either an alternative start or a no codon start outweighs the stringency values

                                # Cull the top hit if it doesn't meet our minimum length requirement anymore, or add it to the temporary list of ORF hits from this nucleotide sequence
                                if len(topHit) < minProLen:                                             # Culling is necessary since we will have shortened the sequence somewhat unless we accepted the topHit as being a no codon start ORF. Note that we will here consider a stop codon in the length of the protein, such that a protein with 99AAs and a stop will pass a minimum 100AA length test. I think this is fair since not all regions here have a stop codon which allows weight to be added to these cases, especially since a stop codon is still conserved as part of an ORF.
                                        doNothing = ''
                                elif maxProLen!= 0 and len(topHit) > maxProLen:                         # Culling should not be necessary here in almost all scenarios, but who knows?
                                        doNothing = ''
                                elif topHit == mPro:                                                    # These temp lists will be populated with potential ORFs from a single nucleotide sequence before being processed in the next major chunk of code starting with 'if len(tempMList + tempAltList + tempNoneList) >= 1:'
                                        if sequenceType.lower() == 'prot':
                                                tempMProt.append(topHit)
                                        elif sequenceType.lower() == 'nucl':
                                                newStartPosition = acceptedPro.find(mPro)
                                                tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                        else:
                                                tempMProt.append(topHit)
                                                newStartPosition = acceptedPro.find(mPro)
                                                tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                elif topHit == altPro:
                                        if sequenceType.lower() == 'prot':
                                                tempAltProt.append(topHit)
                                        elif sequenceType.lower() == 'nucl':
                                                newStartPosition = acceptedPro.find(altPro[1:]) - 1     # - 1 since we're looking at the second character in our altPro (just in case we're replacing with M)
                                                tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                        else:
                                                tempAltProt.append(topHit)
                                                newStartPosition = acceptedPro.find(altPro[1:]) - 1
                                                tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                elif topHit == nonePro:
                                        if sequenceType.lower() == 'prot':
                                                tempNoneProt.append(topHit)
                                        elif sequenceType.lower() == 'nucl':
                                                tempNoneNucl.append(str(nucSeqOfProt))
                                        else:
                                                tempNoneProt.append(topHit)
                                                tempNoneNucl.append(str(nucSeqOfProt))

        # Sort our top hits from each inter-stop codon fragment by size and category (i.e. mPro or altPro?) and select the top X hits
        if len(tempMProt + tempAltProt + tempNoneProt) >= 1 or len(tempMNucl + tempAltNucl + tempNoneNucl) >= 1:
                # Append '-' entries to lists which have less entries than we want to pull to allow the below 'for' loops to run without exceptions
                        # Prot list     [If we are only looking at nucleotides, then prot lists will be populated with hyphens which, realistically, won't impact memory consumption
                for i in range(0, hitsToPull-len(tempMProt)):
                        tempMProt.append('-')
                for i in range(0, hitsToPull-len(tempAltProt)):
                        tempAltProt.append('-')
                for i in range(0, hitsToPull-len(tempNoneProt)):
                        tempNoneProt.append('-')
                        # Nucl list
                for i in range(0, hitsToPull-len(tempMNucl)):
                        tempMNucl.append('-')
                for i in range(0, hitsToPull-len(tempAltNucl)):
                        tempAltNucl.append('-')
                for i in range(0, hitsToPull-len(tempNoneNucl)):
                        tempNoneNucl.append('-')
                # Sort the lists by size (largest on the bottom to allow the .pop() method to remove a hit when accepted) [as above, the prot or nucl variants might just be lists of hyphens. Running the sort twice shouldn't realistically impact time in that case.
                        # Prot
                tempSortedMProt = sorted(tempMProt, key=len)
                tempSortedAltProt = sorted(tempAltProt, key=len)
                tempSortedNoneProt = sorted(tempNoneProt, key=len)
                        # Nucl
                tempSortedMNucl = sorted(tempMNucl, key=len)
                tempSortedAltNucl = sorted(tempAltNucl, key=len)
                tempSortedNoneNucl = sorted(tempNoneNucl, key=len)
                # Run a final size comparison to choose the best ORF(s). We need to split this into two separate statements since the stringency values need to be *3 for nucls
                if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):
                                if len(tempSortedNoneProt[-1]) > len(tempSortedAltProt[-1]) + noCodonStringency and len(tempSortedNoneProt[-1]) > len(tempSortedMProt[-1]) + noCodonStringency:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                        tempOverallProt.append(tempSortedNoneProt[-1])
                                        tempSortedNoneProt.pop()
                                elif len(tempSortedAltProt[-1]) > len(tempSortedMProt[-1]) + altCodonStringency:
                                        tempOverallProt.append(tempSortedAltProt[-1])
                                        tempSortedAltProt.pop()
                                else:
                                        tempOverallProt.append(tempSortedMProt[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                        tempSortedMProt.pop()
                if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):
                                if len(tempSortedNoneNucl[-1]) > len(tempSortedAltNucl[-1]) + noCodonStringency*3 and len(tempSortedNoneNucl[-1]) > len(tempSortedMNucl[-1]) + noCodonStringency*3:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                        tempOverallNucl.append(tempSortedNoneNucl[-1])
                                        tempSortedNoneNucl.pop()
                                elif len(tempSortedAltNucl[-1]) > len(tempSortedMNucl[-1]) + altCodonStringency*3:
                                        tempOverallNucl.append(tempSortedAltNucl[-1])
                                        tempSortedAltNucl.pop()
                                else:
                                        tempOverallNucl.append(tempSortedMNucl[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                        tempSortedMNucl.pop()
                # Format and produce the output of this script
                tempOutputProt = []
                tempOutputNucl = []
                if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):                          
                                if tempOverallProt[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                        break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                                tempOutputProt.append('>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallProt[i])
                        if len(tempOutputProt) > 0:    # Need this check for sequences that don't get any hits
                                outputProt.append('\n'.join(tempOutputProt))

                if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                        for i in range(0, hitsToPull):                          
                                if tempOverallNucl[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                        break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                                tempOutputNucl.append('>' + record.id + '_ORF' + str(i+1) + '\n' + tempOverallNucl[i])
                        if len(tempOutputNucl) > 0:
                                outputNucl.append('\n'.join(tempOutputNucl))
        else:
                doNothing = ''                # We don't need to do anything if no hits were found that pass the minimum length threshold

        ongoingCount += 1
        
        # Save backup if ongoingCount == 10,000. There are two sections here, the first will create the file on the first loop, the second will add to the file on subsequent loops
        if sequenceType.lower() == 'prot':
                if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                        with open(outputFileName, 'w') as output:
                                output.write('\n'.join(outputProt))
                        outputProt = []
                elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                        with open(outputFileName, 'a') as output:
                                output.write('\n')
                                output.write('\n'.join(outputProt))
                        outputProt = []
        elif sequenceType.lower() == 'nucl':
                if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                        with open(outputFileName, 'w') as output:
                                output.write('\n'.join(outputNucl))
                        outputNucl = []
                elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                        with open(outputFileName, 'a') as output:
                                output.write('\n')
                                output.write('\n'.join(outputNucl))
                        outputNucl = []
        else:
                if ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + protOutName) == False:       # Doesn't matter if we check for protOutName or nuclOutName. Theoretically, if one of the files is deleted while the program is running it could be problematic, but I mean, what can I really do about that without child-proofing the script excessively?
                        with open(protOutName, 'w') as protFile, open(nuclOutName, 'w') as nuclFile:
                                protFile.write('\n'.join(outputProt))
                                nuclFile.write('\n'.join(outputNucl))
                        outputProt = []
                        outputNucl = []
                elif ongoingCount%10000 == 0 and os.path.isfile(os.getcwd() + '\\' + protOutName) == True:
                        with open(protOutName, 'a') as protFile, open(nuclOutName, 'a') as nuclFile:
                                protFile.write('\n')
                                nuclFile.write('\n')
                                protFile.write('\n'.join(outputProt))
                                nuclFile.write('\n'.join(outputNucl))
                        outputProt = []
                        outputNucl = []

# Dump the last few results after the script has finished, or create the output if there were less than 10,000 sequences
if sequenceType.lower() == 'prot':
        if os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                with open(outputFileName, 'w') as output:
                        output.write('\n'.join(outputProt))
        elif os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                with open(outputFileName, 'a') as output:
                        output.write('\n')
                        output.write('\n'.join(outputProt))
elif sequenceType.lower() == 'nucl':
        if os.path.isfile(os.getcwd() + '\\' + outputFileName) == False:
                with open(outputFileName, 'w') as output:
                        output.write('\n'.join(outputNucl))
        elif os.path.isfile(os.getcwd() + '\\' + outputFileName) == True:
                with open(outputFileName, 'a') as output:
                        output.write('\n')
                        output.write('\n'.join(outputNucl))
else:
        if os.path.isfile(os.getcwd() + '\\' + protOutName) == False:       # Doesn't matter if we check for protOutName or nuclOutName. Theoretically, if one of the files is deleted while the program is running it could be problematic, but I mean, what can I really do about that without child-proofing the script excessively?
                with open(protOutName, 'w') as protFile, open(nuclOutName, 'w') as nuclFile:
                        protFile.write('\n'.join(outputProt))
                        nuclFile.write('\n'.join(outputNucl))
        elif os.path.isfile(os.getcwd() + '\\' + protOutName) == True:
                with open(protOutName, 'a') as protFile, open(nuclOutName, 'a') as nuclFile:
                        protFile.write('\n')
                        nuclFile.write('\n')
                        protFile.write('\n'.join(outputProt))
                        nuclFile.write('\n'.join(outputNucl))
exit()







