#!/usr/bin/env python

"""
# pA-TSD_clearner.py
# usage: python %prog infile.txt
# python2.7
"""

import os,sys,re

f_path=sys.argv[1]
basename,ext=os.path.splitext(f_path)
outfile_path=basename+'_pA_clean'+ext

# 3x-long AT-rich simple repeats defined by RepeatMasker
rep='acacac|atatat|agagag|aataataat|ataataata|aacaacaac|acaacaaca|aagaagaag|agaagaaga|agatagatagat|atagatagatag|aaacaaacaaac|aacaaacaaaca|acaaacaaacaa|aaataaataaat|aataaataaata|ataaataaataa|aaagaaagaaag|aagaaagaaaga|agaaagaaagaa|acatacatacat|atacatacatac|acagacagacag|agacagacagac|aaaacaaaacaaaac|aaacaaaacaaaaca|aacaaaacaaaacaa|acaaaacaaaacaaa|aaaataaaataaaat|aaataaaataaaata|aataaaataaaataa|ataaaataaaataaa|aaaagaaaagaaaag|aaagaaaagaaaaga|aagaaaagaaaagaa|agaaaagaaaagaaa|acaatacaatacaat|aatacaatacaatac|atacaatacaataca|agaatagaatagaat|aatagaatagaatag|atagaatagaataga|atacatatacatatacat|acatatacatatacatat|atatacatatacatatac|aaaaacaaaaacaaaaac|aaaacaaaaacaaaaaca|aaacaaaaacaaaaacaa|aacaaaaacaaaaacaaa|acaaaaacaaaaacaaaa|aaaaagaaaaagaaaaag|aaaagaaaaagaaaaaga|aaagaaaaagaaaaagaa|aagaaaaagaaaaagaaa|agaaaaagaaaaagaaaa|aaaaataaaaataaaaat|aaaataaaaataaaaata|aaataaaaataaaaataa|aataaaaataaaaataaa|ataaaaataaaaataaaa|agacagagacagagacag|acagagacagagacagag|agagacagagacagagac|tgtgtg|tctctc|tatata|attattatt|tattattat|gttgttgtt|tgttgttgt|cttcttctt|tcttcttct|ctatctatctat|atctatctatct|ttgtttgtttgt|tgtttgtttgtt|gtttgtttgttt|ttatttatttat|tatttatttatt|atttatttattt|ttctttctttct|tctttctttctt|ctttctttcttt|gtatgtatgtat|atgtatgtatgt|ctgtctgtctgt|gtctgtctgtct|tttgttttgttttgt|ttgttttgttttgtt|tgttttgttttgttt|gttttgttttgtttt|tttattttattttat|ttattttattttatt|tattttattttattt|attttattttatttt|tttcttttcttttct|ttcttttcttttctt|tcttttcttttcttt|cttttcttttctttt|tgtattgtattgtat|gtattgtattgtatt|attgtattgtattgt|tctattctattctat|ctattctattctatt|attctattctattct|atgtatatgtatatgtat|gtatatgtatatgtatat|atatgtatatgtatatgt|ttttgtttttgtttttgt|tttgtttttgtttttgtt|ttgtttttgtttttgttt|tgtttttgtttttgtttt|gtttttgtttttgttttt|ttttctttttctttttct|tttctttttctttttctt|ttctttttctttttcttt|tctttttctttttctttt|ctttttctttttcttttt|ttttatttttatttttat|tttatttttatttttatt|ttatttttatttttattt|tatttttatttttatttt|atttttatttttattttt|ctgtctctgtctctgtct|gtctctgtctctgtctct|ctctgtctctgtctctgt'

outfile=open(outfile_path, 'w')

with open(f_path, 'r') as infile:
    next(infile)
    lines=['polyA.or.polyT\tTSD_query\tTSD_search\tline_number\n']
    pA=''
    n=0
    for line in infile:
        if re.search(r'%s' % rep, line) is None:
            hit=re.split(r'\t|,', line)
            if (len(hit[3]) == 7 and hit[3].count('A') < 5):
                if hit[1] == pA:
                    n=n+1
                    lines.append(','.join(hit).replace(',>', '\t>'))
                else:
                    pA=hit[1]
                    if n < 4:
                        outfile.write(''.join(lines))
                    lines=[]
                    lines.append(','.join(hit).replace(',>', '\t>'))
                    n=1
            elif (len(hit[3]) == 9 and hit[3].count('A') < 6):
                if hit[1] == pA:
                    n=n+1
                    lines.append(','.join(hit).replace(',>', '\t>'))
                else:
                    pA=hit[1]
                    if n < 4:
                        outfile.write(''.join(lines))
                    lines=[]
                    lines.append(','.join(hit).replace(',>', '\t>'))
                    n=1
            elif (len(hit[3]) == 11 and hit[3].count('A') < 7):
                if hit[1] == pA:
                    n=n+1
                    lines.append(','.join(hit).replace(',>', '\t>'))
                else:
                    pA=hit[1]
                    if n < 4:
                        outfile.write(''.join(lines))
                    lines=[]
                    lines.append(','.join(hit).replace(',>', '\t>'))
                    n=1
    if n < 4:
        outfile.write(''.join(lines))

outfile.close()
