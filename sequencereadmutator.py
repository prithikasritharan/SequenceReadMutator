'''
21/05/2018
Python script to generate mutations at different positions in DNA sequence reads
'''

import sys
import argparse
import random
from itertools import cycle, islice, chain

#Instantiate the parser
parser = argparse.ArgumentParser(description='Generate mutations at different positions in DNA sequence reads')

#Required positional argument
parser.add_argument('seq_file', type=str,
                    help='A required argument stating the sequence reads file')
parser.add_argument('out_file', type=str,
                    help='Output file to write mutated sequences')

#Parse arguments
args = parser.parse_args()



#function that takes the length of sequence as input and uses the random function to randomly choose A, C, G or T for the length of the sequence
def generate_string(N, alphabet='ACGT'):
    return ''.join([random.choice(alphabet) for i in xrange(N)])


#function that calculates the maximum number of mutations to be induced within a read and cycles this list so that it is the same length as the number of reads to be mutated
def no_of_mutations(a, l):
    #Calculate the maximum number of mutations to be induced by multiplying i by then read length
    mut_length = int(a*read_length)
    #create an array from 1 to mut_length (the max no. of mutations)
    mutations = range(1, (mut_length+1))
    #print 'The maximim number of mutations to be induced is: ' + str(mutations[-1])
    
    #cycles through the no. of mutations in a list the same length as the no. of reads 
    cycled = cycle(mutations)
    sliced_mut = islice(cycled, None, l)
    return sliced_mut


#function that mutates the starting bases of reads by taking a list of reads as input 
def mutate_start(first_mut):
    length = len(first_mut) #takes length of list containing the reads
    new_reads = [] #initialise list to add the mutated reads
    #calls the no_of_mutations function which takes the value of m (max no of mut) and no of reads and returns a list containing the no of mutations to be induced 
    mut_array = list(no_of_mutations(m, length))
    for i in range(length):
        #loops through the reads and converts each read into a list
        read = list(first_mut[i])
        #generates random mutation using the generate_string function based on the no of mutations at the corresponding pos. in the mut_array list 
        var = list(generate_string(mut_array[i]))
        #compares the bases of the mutation with the read base to ensure they are not the same before altering the base sequence of the read at that position to the mutated base
        for x, base in enumerate(var):
            while base == read[x]:
                base = generate_string(1)
                if base != read[x]:
                    break
            read[x] = base
        #converts bases of the read back into a string and appends read to the new_reads list which is returned
        fin_read = ''.join(read)
        new_reads.append(fin_read[:])
    return mut_array, new_reads #returns a list containing the no of mutations in each read and the mutated read sequences


#function that mutates the last few bases of reads by taking a list of reads as input
def mutate_end(end_mut):
    length = len(end_mut) #takes length of list containing the reads
    new_reads = [] #initialise list to add the mutated reads
    mut_array = list(no_of_mutations(m, length))
    for i in range(length):
        #loops through the reads, reverses the read sequence and convert it into a list
        read = list(reversed(end_mut[i]))
        var = list(generate_string(mut_array[i]))
        #compares mutated base with read base to ensure they are not the same before changing the value of the read base to the mutated base 
        for x, base in enumerate(var):
            while base == read[x]:
                base = generate_string(1)
                if base != read[x]:
                    break
            read[x] = base
        rev_read = reversed(read) #reverses the read sequence back to the initial order (minus the mutated bases)
        fin_read = ''.join(rev_read)
        new_reads.append(fin_read[:])
    return mut_array, new_reads #returns a list containing the no of mutations in each read and the mutated read sequences  


#function that mutates a few bases at a random position in a read by taking a list of reads as input
def mutate_clump(third_mut):
    length = len(third_mut)
    new_reads = []
    mut_array = list(no_of_mutations(m, length))
    for i in range(length):
        read = list(third_mut[i])
        var = list(generate_string(mut_array[i]))
        #generates a random int between 0 and the length of the read minus the no. of mutations to be induced 
        rand_pos = random.randint(0,(len(read)-len(var)))
        #j is the intial starting position of the base to be mutated in the read
        j = rand_pos
        #compares mutated base with the read base starting at position j to ensure they are not the same before changing the value of the read base to the mutated base
        for base in var:
            while base == read[j]:
                base = generate_string(1)
                if base != read[j]:
                    break
            read[j] = base
            j += 1
        fin_read = ''.join(read)
        new_reads.append(fin_read[:])
    return mut_array, new_reads


#function that mutates bases at many random positions in a read by taking a list of reads as input
def mutate_random(fourth_mut):
    length = len(fourth_mut)
    new_reads = []
    mut_array = list(no_of_mutations(m, length))
    for i in range(length):
        read = list(fourth_mut[i])
        var = list(generate_string(mut_array[i]))
        pos = random.sample(range(len(read)), mut_array[i]) #generates a list of random positions to be mutated so that there aren't any duplicate position values
        for index, j in enumerate(pos):
            while var[index] == read[j]: #compares the sequence at the random position in read to base
                var[index] = generate_string(1)
                if var[index] != read[j]:
                    break
            read[j] = var[index]            
        fin_read = ''.join(read)
        new_reads.append(fin_read[:])
    return mut_array, new_reads


#function that returns a list with the different types of mutations across the reads to write out
def list_of_mutations():
    mut_type = ["No mutation"] * len(no_mut)
    mut_type.extend(["Start mutation"] * len(start_mut))
    mut_type.extend(["End mutation"] * len(end_mut))
    mut_type.extend(["Clump mutation"] * len(clump_mut))
    mut_type.extend(["Random mutation"] * len(rand_mut))
    return mut_type  


#read in the sequence reads file into a list
seq_file = open(args.seq_file, "r")
reads = seq_file.read().splitlines()
seq_file.close()

#get the number of reads in the list
read_num = len(reads)
print 'The number of reads in the file is: ' + str(read_num)

#get read length
read_length = len(reads[0])
print 'The length of each read is: ' + str(read_length) + ' bases'

#set a value for m which will determine the maximum number of mutations
m = 0.3

#percentage of reads with different types of mutations
no_mut = reads[0:int(0.1*read_num)] #slices the reads that aren't going to be mutated
start = int(0.2*read_num) #no. of reads with mutation at start
end_pos = len(no_mut) + start #starting position of end mutations
start_mut = reads[len(no_mut):end_pos] #slices the no. of reads that are going to have start mutations
end = int(0.2*read_num) #reads with mutation at end
clump_pos = end_pos + end #starting position of clumped mutations
end_mut = reads[end_pos:clump_pos] #slices the no. of reads that are going to have end mutations
clump = int(0.2*read_num) #reads with mutations together starting at a random position
rand_pos = clump_pos + clump #starting position of random mutations
clump_mut = reads[clump_pos:rand_pos] #slices the no. of reads that are going to have mutations that are clumped together
rand = int(0.3*read_num) #reads with mutations at random positions
rand_mut = reads[rand_pos:] #slices the no. of reads that are going to have random mutations


no0 = [0] * len(no_mut)
#mutates the sequences at the start of the read by calling the mutate_start function on a set of reads
no1, mut1 = mutate_start(start_mut)
#mutates the sequences at the end of the read by calling the mutate_end function on a set of reads
no2, mut2 = mutate_end(end_mut)
#mutates a few bases together at a random position in the read by calling the mutate_clump function on a set of reads
no3, mut3 = mutate_clump(clump_mut)
#mutates bases at random positions throughout the read by calling the mutate_random function on a set of reads
no4, mut4 = mutate_random(rand_mut)

#join the lists containing the no.of mutations across all the different functions
final_nom = chain(no0, no1, no2, no3, no4)
#join the lists containing the mutated reads across all the different functions
final_read = chain(no_mut, mut1, mut2, mut3, mut4)

#obtain a list of the type of mutation across the reads
mutation_type = list_of_mutations()
 

#output file to write mutated sequences and second file that contains the sequences and the mutations within that read
with open(args.out_file, "w") as f, open("{0}.withmutations".format(args.out_file), "w") as nf: 
    #Get the first set of reads without any mutations and write it to output file
        for line, mut_line, type_line in zip(final_nom, final_read, mutation_type):
            f.write("%s\n" % mut_line)
            nf.write("%s\t%s\t%s\n" % (mut_line, line, type_line))


