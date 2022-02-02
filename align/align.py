# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0
        self.score=[]

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        Needleman-Wunsch Algorithm 
         performs global sequence alignment of seqA and seqB
        
        Parameters:
            seqA: Sequence A to aling
            seqB: Sequence B to align 
        Returns:
        tuple with the following format:
        (alignment score, seqA alignment, seqB alignment)
        """
        
        ## I used this values to intialize the matrices
        Alen=len(seqA) + 1
        Blen=len(seqB)+1
        SeqA = [char for char in seqA]
        SeqB = [char for char in seqB]

        # Initialize 6 matrix private attributes for use in alignment
        # create matrices for alignment scores and gaps
        self._align_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        ## This matrix in the position 0,0 is always going to be 0
        self._align_matrix[0,0]=0
        ## Now matrix for gap in A 
        self._gapA_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        ## GAP A starts with a gap open  
        self._gapA_matrix[0,0]=self.gap_open
        #Then we add a gap extension 
        for i in range(Alen-1):
            i=i+1
            self._gapA_matrix[i,0]=self._gapA_matrix[i-1,0]+self.gap_extend
        
        # This matrix also starts with a gap openning and then adding the gap extend   
        self._gapB_matrix = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._gapB_matrix[0,0]=self.gap_open
        for i in range(Blen-1):
            i=i+1
            self._gapB_matrix[0,i]=self._gapB_matrix[0,i-1]+self.gap_extend

        # create matrices for pointers used in backtrace procedure
        self._back = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_A = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        self._back_B = np.ones((len(seqA) + 1, len(seqB) + 1)) * -np.inf
        ## The  pointers  for GAP A for the first column are always going to point up
        for i in range(Alen):
            self._back_A[i,0]=1
        ## for B the pointer the first row are always going to go left
        for i in range(Blen):
            self._back_B[0,i]=2
        
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB

        # Starting the global sequence alignment 
        for i in range(Alen-1):
            i=i+1 ##skip the first colum since is already filled for all the matrices
            for j in range(Blen-1):
                j=j+1 #skip the first row since is already filled
                
                ## Matrix align
                scoreM=self._align_matrix[i-1,j-1] ## compute score from aligned matrix
                scoreX=self._gapA_matrix[i-1,j-1] ## compute score from A matrix
                scoreY=self._gapB_matrix[i-1,j-1] ## comoute score from B matrix
        
                max_value = max(scoreM,scoreX,scoreY) ##choose the max value
                All_scores = [scoreM,scoreX,scoreY]
                ## Add to the pointers where the max value came from
                self._back[i,j]=All_scores.index(max_value)
                ##Add the max value plus the matching score in that position
                self._align_matrix[i,j] = self.sub_dict[SeqA[i-1],SeqB[j-1]]+ max_value
    
                ## Matrix A
                scoreM=(self.gap_open + self.gap_extend + self._align_matrix[i-1,j]) #score for align matrix
                scoreX=( self.gap_extend + self._gapA_matrix[i-1,j]) ## score for matrix gap A
                scoreY=(self.gap_open +  self.gap_extend + self._gapB_matrix[i-1,j]) ## score for matrix gap B
        
                max_value=max(scoreM,scoreX,scoreY) ##choose the max value to add to the matrix
                All_scores = [scoreM,scoreX,scoreY]
                self._back_A[i,j]=All_scores.index(max_value) ## add the index of the max value to the pointer
                self._gapA_matrix[i,j]=max_value ## add the max to the gapA matrix
    
                #Matrix B
                scoreM=(self.gap_open + self.gap_extend + self._align_matrix[i,j-1]) # compute score for matrix align
                scoreX=(self.gap_open + self.gap_extend + self._gapA_matrix[i,j-1]) ## compute score for matrix A
                scoreY=(self.gap_extend + self._gapB_matrix[i,j-1])#score for matrix B
        
                max_value=max(scoreM,scoreX,scoreY)
                All_scores = [scoreM,scoreX,scoreY]
                self._back_B[i,j]=All_scores.index(max_value) ## add the index of the max value to the pointer matrix
                self._gapB_matrix[i,j] = max_value # add the max value to the B matrix


        return self._backtrace()

    def M_function(self,row,col):
        """

        Parameters
        ----------
        row : int
            row we are currently traversing.
        col : int
            colum we are currently traversing.

        Returns
        -------
        If reached the end of the matrix return score
        If not then keep looking at the appropiate pointer
        
        """
        if col==0 and row==0:
            return(self.score)
        else:
            ##since we are at align matrix this means we match the value for both strings
            self.score.append(1)
            pointer_value = self._back[row,col]
            #0 = keep looking at this matrix for pointer in the diagonal position
            if pointer_value==0:
                return(self.M_function(row-1,col-1))
            # 1 = look at value gapA matrix and back_gapA (diagonal position)
            if pointer_value==1:
                return(self.X_function(row-1,col-1))
            #2 = look at pointer gapB matrix and back_gapB (diagonal position)
            if pointer_value==2:
                return(self.Y_function(row-1,col-1))         
    
    def X_function(self,row,col):
        """ 
        Parameters
        ----------
        row : int
            row we are currently traversing.
        col : int
            colum we are currently traversing.

        Returns
        -------
        If reached the end of the matrix return score
        If not then keep looking at the appropiate pointer

        """
        if col==0 and row==0:
            return(self.score)
        else:
       ## since we are at the gap A matrix this means we have a gap with regards of sequence B
            self.score.append("GAP IN B")
            pointer_value = self._back_A[row,col]
            # 0=look the _back matrix in the up postion
            if pointer_value==0:
                return(self.M_function(row-1,col))
            # 1=look the _back_gapA matrix in the up postion
            if pointer_value==1:
                return(self.X_function(row-1,col))
            # 2=look the _back_gapB matrix in the up postion
            if pointer_value==2:
                return(self.Y_function(row-1,col))          
    def Y_function(self,row,col):
        """ 
        Parameters
        ----------
        row : int
            row we are currently traversing.
        col : int
            colum we are currently traversing.

        Returns
        -------
        If reached the end of the matrix return score
        If not then keep looking at the appropiate pointer

        """
        if col==0 and row==0:
            return(self.score)
        else:
            ## since we are at the gap B matrix this means we have a gap with regards of sequence A
            self.score.append("GAP IN A")
            pointer_value = self._back_B[row,col]
            # 0=look the _back matrix in the left postion
            if pointer_value==0:
                return(self.M_function(row,col-1))
            # 1=look the _back_gapA matrix in the left
            if pointer_value==1:
                return(self.X_function(row,col-1))
            # 1=look the _back_gapB matrix in the left
            if pointer_value==2:
                return(self.Y_function(row,col-1))  
            
    def _backtrace(self) -> Tuple[float, str, str]:
        """
        Method to recover the path for the sequence alignment
        Returns to align the sequece score and the string for alignment
        It uses the recursive methods: M_function, X_function and Y_function to keep track of the pointers
        
        Attributes:
            score: list
            List of the sequence for alignmet_seqA and alignment_seqB
        
        """
        col = len(self._seqB)
        row = len(self._seqA)
        ## We check the three matrices at the last position and select the maximum score
        max_start = max(self._align_matrix[row,col], self._gapA_matrix[row,col], self._gapB_matrix[row,col])
        
        self.alignment_score = max_start ## here we start traversing and also add the final score
        final_scores = [self._align_matrix[row,col], self._gapA_matrix[row,col], self._gapB_matrix[row,col]]
        ## get the index so we can start at the corret pointer
        max_index=final_scores.index(max_start) 
        
        #If the max value is in postion 0 then start at the pointer for the aling matrix
        if max_index==0:
            self.M_function(row,col)
        #If the max value is in postion 1 then start at the pointer for the gap A matrix
        if max_index==1:
            self.X_funtion(row,col)
        #If the max value is in postion 2 then start at the pointer for the gap B matrix
        if max_index==2:
            self.Y_function(row,col)
        
        ## strings to calculate the last seq alignment
        SA=""
        SB=""
        # Used as a temporary list before producing final alignment  
        SA_Alignment=[]
        SB_Alignment=[]
        
        ## Sepaarte the GAPS in A and GAPS in B
        for i in (self.score):
            if i==1:
                SA_Alignment.append(1)
                SB_Alignment.append(1)
            elif i == "GAP IN A":
                SA_Alignment.append("-")
                SB_Alignment.append(1)
            elif i == 'GAP IN B':
                SB_Alignment.append("-")
                SA_Alignment.append(1)  

        seqA_reverse = self._seqA[::-1] 
        seqB_reverse = self._seqB[::-1] 
        
        ## For each string if find one print that letter, if not add a gap
        count=0
        for i in SA_Alignment:
            if i==1:
                SA+=seqA_reverse[count]
                count+=1
            else:
                SA+="-"
        count=0
        for i in SB_Alignment:
            if i==1:
                SB+=seqB_reverse[count]
                count+=1
            else:
                SB+="-"    
        ## Reverse the string and  use it as a antribute
        self.seqA_align = SA[::-1]
        self.seqB_align = SB[::-1]
        return(self.alignment_score,self.seqA_align,self.seqB_align)
            


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header

