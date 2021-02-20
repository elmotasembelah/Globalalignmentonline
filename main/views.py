from django.http import QueryDict
from django.shortcuts import render, redirect
import numpy as np


# Create your views here.


def homepage(request):
    if request.method == 'GET':
        return render(request, template_name='main/home.html')

    elif request.method == 'POST':
        userinput = QueryDict.copy(request.POST)

        seq1 = userinput['seq1']
        seq2 = userinput['seq2']

        result = dnacheck(seq1, seq2)

        return render(request, template_name='main/result.html', context={'resultdict': result})


def resultpage(request):
        return render(request, template_name='main/result.html')



def dnacheck(seq1, seq2):
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    nucleotides = 'ACGT'
    check = False

    for char in seq1:
        if nucleotides.find(char) == -1:
            check = True
            break

    if not check:
        for char in seq2:
            if nucleotides.find(char) == -1:
                check = True

    if check:
        return proteinsequencecheck(seq1, seq2)
    else:
        return algorithm(seq1, seq2)


def proteinsequencecheck(seq1, seq2):
    resultdict = {}
    seq1 = seq1.upper()
    seq2 = seq2.upper()

    components = 'ARNDCQEGHILKMFPSTWYV'
    check = False

    for char in seq1:
        if components.find(char) == -1:
            check = True
            break

    if not check:
        for char in seq2:
            if components.find(char) == -1:
                check = True

    if check:
        resultdict['alig1'] = 'One or both the sequences that you entered were neither DNA or a protein sequence, please try again after reviewing the sequences'
        return resultdict

    else:
        return algorithm(seq1, seq2)


def algorithm(seq1, seq2):
    resultdict = {}

    main_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    match_checker_matrix = np.zeros((len(seq1), len(seq2)))

    match_reward = 1
    mismatch_penalty = -1
    gap_penalty = -2

    # Fill the match checker matrix accrording to match or mismatch
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                match_checker_matrix[i][j] = match_reward
            else:
                match_checker_matrix[i][j] = mismatch_penalty


    for i in range(len(seq1) + 1):
        main_matrix[i][0] = i * gap_penalty
    for j in range(len(seq2) + 1):
        main_matrix[0][j] = j * gap_penalty

    # STEP 2 : Matrix Filling
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            main_matrix[i][j] = max(main_matrix[i - 1][j - 1] + match_checker_matrix[i - 1][j - 1],
                                    main_matrix[i - 1][j] + gap_penalty,
                                    main_matrix[i][j - 1] + gap_penalty)



    aligned_1 = ""
    aligned_2 = ""

    ti = len(seq1)
    tj = len(seq2)

    while ti > 0 and tj > 0:

        if ti > 0 and tj > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj - 1] + match_checker_matrix[ti - 1][tj - 1]:
            aligned_1 = seq1[ti - 1] + aligned_1
            aligned_2 = seq2[tj - 1] + aligned_2

            ti = ti - 1
            tj = tj - 1

        elif ti > 0 and main_matrix[ti][tj] == main_matrix[ti - 1][tj] + gap_penalty:
            aligned_1 = seq1[ti - 1] + aligned_1
            aligned_2 = "-" + aligned_2

            ti = ti - 1
        else:
            aligned_1 = "-" + aligned_1
            aligned_2 = seq2[tj - 1] + aligned_2

            tj = tj - 1

    resultdict['alig1'] = aligned_1
    resultdict['alig2'] = aligned_2



    return resultdict
