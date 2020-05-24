#!/usr/bin/env python
# sends a combined email for all the fasttests run with different versions of matlab

# Reference
# https://stackoverflow.com/questions/3362600/how-to-send-email-attachments
import pdb
import sys
import os
import glob
import smtplib
import datetime
import subprocess
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate
#from email.headerregistry import Address

# change these email addresses if needed
to_emails = ['Sou-Cheng Terrya Choi <schoi32@iit.edu>', 'Sou-Cheng Terrya Choi <sou.cheng.terrya.choi@gmail.com>']
cc_emails = ['Jagadeeswaran <jrathin1@iit.edu>', 'Kan Zhang <kzhang23@hawk.iit.edu>']  # 'GAIL Dev <gail_dev@googlegroups.com>',

# to_emails = ['Jagadeeswaran <jrathin1@iit.edu>' ]
# cc_emails = ['Jagadeeswaran <jagadeesr@gmail.com>' ]

# dont chnage this
from_email = 'GAIL Project <gail@local.iit.edu>'

def send_mail(send_from, send_to, send_cc, subject, text, files=None,
              server="127.0.0.1"):
    assert isinstance(send_to, list)
    assert isinstance(send_cc, list)

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = COMMASPACE.join(send_to)
    msg['Cc'] = COMMASPACE.join(send_cc)
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject

    msg.attach(MIMEText(text, "plain", "utf-8"))

    for f in files or []:
        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=basename(f)
            )
            part['Content-Disposition'] = 'attachment; filename="%s.txt"' % basename(f)
            msg.attach(part)

    smtp = smtplib.SMTP(server)
    smtp.sendmail(send_from, send_to+send_cc, msg.as_string())
    smtp.close()

datetime_today = datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
subject = "Daily test results: {0}".format(datetime_today)
text = 'Please review the output log for more details'
date_today = datetime.datetime.now().strftime("%Y-%m-%d")

# path where all the error greps stored
path_automatic_tests = '/home/gail/GAIL_tests/repo/gail-development/GAIL_Matlab/Tests/Developers_only/automatic_tests/'
Test_results_correct = 'Test_Results_correct.txt'
Test_results_file_prefix = 'Test_Results-'

def prcoess_test_results(path, correct_file, Test_results_file_prefix):
    # compare the result file (Test_ResultsR2016a.tx) with the reference-correct file Test_Results_correct.txt
    search_pattern = path + Test_results_file_prefix + 'R*.txt'
    test_result_files = glob.glob(search_pattern)
    compare_result = {}
    for file_t in test_result_files:
        fn = basename(file_t)
        fn = fn.strip('.txt')
        matlabVer = fn.replace(Test_results_file_prefix, '')
        with open(path_automatic_tests+correct_file, 'r') as ref_file:
            with open(file_t, 'r') as test_file:
                same = set(test_file).issubset(ref_file)

        if same is False:
            compare_result.update({matlabVer:('Fail', file_t)})
        else:
            compare_result.update({matlabVer:('Pass', None)})

    return compare_result

# log_output_file="/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-$(date +\%Y-\%m-\%d).out"

if __name__ == '__main__':
    # pdb.set_trace()

    log_output_file = sys.argv[1]
    compare_result = prcoess_test_results(path_automatic_tests, Test_results_correct, Test_results_file_prefix)

    mail_body = 'We tested with the following Matlab versions \n\n'
    mail_body += 'Matlab : Result\n'
    mail_body += '---------------------------\n'
    output_files = [log_output_file]
    tests_passed = True
    error_text = ''

    # check if pass/fail then create message to pring and collect error greps
    for m_ver, res in sorted(compare_result.items()):
        if res[0] is 'Fail':
            # If the test failed then attach the error grep
            with open(res[1], 'r') as test_file:
                error_text += m_ver + ': \n' + test_file.read() + '\n'
            tests_passed = False
        mail_body += m_ver + ' : ' + res[0] + '\n'

    mail_body += '\n' + error_text + '\n' + text
    mail_body += '\n'
    def mega_cmd_login():
        MEGASYNC_UNAME='jrathin1@iit.edu'
        MEGASYNC_PWD='TestsForGail'
        proc = subprocess.Popen(['/home/gail/usr/bin/mega-login', MEGASYNC_UNAME, MEGASYNC_PWD],stdout=subprocess.PIPE)
        if True:
            line = proc.stdout.readline()
            if line != '':
                # got output: parse it
                print("Result:", line)

                
    def get_file_link(fn):
        result = ''
        proc = subprocess.Popen(['/home/gail/usr/bin/mega-export','-a', fn],stdout=subprocess.PIPE)
        if True:
            line = proc.stdout.readline().decode('utf-8')
            if line != '':
                # got output: parse it
                print("Result:", line.rstrip())
                if 'Exported' in line:
                    link = line.rstrip('\n').split(': ')
                    result = line
        
        return result
        
    # mega_cmd_login()
    for fn in output_files:
        filepath = '/PBS_jobs/pbs_reports/' + os.path.split(fn)[-1]
        exp_result = get_file_link(filepath)
        mail_body += "{0} : {1} ".format(os.path.split(fn)[-1], exp_result.split(': ')[1])
    
    mail_body += '\n'
    mail_body += 'Archived reports can be found here \n'
    test_report_dirs = {'PBS logs':'/PBS_jobs/pbs_reports', 'Fast tests': '/test_reports', 'Workouts': '/workout_reports'}
    for tn, trd in test_report_dirs.items():
        exp_result = get_file_link(trd)
        temp = "{0} : {1}".format(tn, exp_result.split(': ')[1])
        mail_body += temp

    subject = "Daily test results {0} : {1}".format( ('OK' if tests_passed else 'Wrong'), datetime_today)
    res = send_mail(from_email, to_emails, cc_emails, subject, mail_body, [])
    print('done')
