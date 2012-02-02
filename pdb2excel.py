__author__ = 'xiaohan'

import xlwt
import os
import re
def toexcel(par_path,o_file):
    wbk =xlwt.Workbook()
    sheet = wbk.add_sheet('sheet1')
    for row_c,complex_id in enumerate(os.listdir(par_path)):
        pattern_file = os.path.join(par_path,complex_id,'pattern.dat')
        with open(pattern_file,'r') as f:
             fp = f.readline().split(':')[1]
             fp = re.sub(r'1\.0','1',fp)

             sheet.write(row_c,0,complex_id)
             sheet.write(row_c,1,fp)

    wbk.save(o_file)

def totext(par_path,o_file):
    with open(o_file,'w') as o_file:
        for row_c,complex_id in enumerate(os.listdir(par_path)):
            pattern_file = os.path.join(par_path,complex_id,'pattern.dat')
            with open(pattern_file,'r') as f:
                 fp = f.readline().split(':')[1].strip()
                 fp = re.sub(r'1\.0','1',fp)
                 o_file.write(complex_id+'\t')
                 for start_id in xrange(len(fp)/18):
                     o_file.write('%s    ' %' '.join(fp[start_id * 18:(start_id + 1) * 18]))
                 o_file.write('\n')
if __name__=='__main__':
    totext('test_output','test_output.dat')

