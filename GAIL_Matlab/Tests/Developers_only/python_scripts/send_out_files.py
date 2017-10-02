#!/usr/bin/env python

# Reference
# https://stackoverflow.com/questions/3362600/how-to-send-email-attachments
import pdb
import os
import smtplib
import datetime
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate
#from email.headerregistry import Address
import zipfile

def send_mail(send_from, send_to, subject, text, files=None,
              server="127.0.0.1"):
    #pdb.set_trace()
    assert isinstance(send_to, list)

    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = COMMASPACE.join(send_to)
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject

    msg.attach(MIMEText(text))

    for f in files or []:
        with open(f, "rb") as fil:
            part = MIMEApplication(
                fil.read(),
                Name=basename(f)
            )
            part['Content-Disposition'] = 'attachment; filename="%s"' % basename(f)
            msg.attach(part)


    smtp = smtplib.SMTP(server)
    smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.close()

def zipdir(path, ziph, pattern):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            if not pattern in file:
                continue
            ziph.write(os.path.join(root, file))

from_email = 'GAIL Project <gail@local.iit.edu>'
to_emails = ['Kan Zhang <kzhang23@hawk.iit.edu>', 'Jagadeeswaran <jrathin1@iit.edu>', 'Sou-Cheng Terrya Choi <schoi32@iit.edu>' ]

datetime_today = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
subject = "Archive of test Output files: {0}".format(datetime_today)
text = 'Please find the attached archive of test output files for the last 30 days'
date_today = datetime.datetime.now().strftime("%Y-%m-%d")
output_files = ['/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-{0}.out'.format(date_today),
                '/home/gail/GAIL_tests/PBS_jobs/pbs_reports/gail_daily_tests-{0}.err'.format(date_today)]
zipped_output_file = '/home/gail/GAIL_tests/PBS_jobs/output_files_all-{0}.out'.format(date_today)
#xx = send_mail('GAIL Project <gail@local.iit.edu>', ['Jagadeeswaran <jrathin1@iit.edu>'], 'Output files', 'hello, testing')
if __name__ == '__main__':
    zipf = zipfile.ZipFile(zipped_output_file, 'w', zipfile.ZIP_DEFLATED)
    zipdir('/home/gail/GAIL_tests/PBS_jobs/pbs_reports/', zipf, 'gail_' )
    zipf.close()
    res = send_mail(from_email, to_emails, subject, text, [zipped_output_file])
    os.remove(zipped_output_file)
    print('done')
