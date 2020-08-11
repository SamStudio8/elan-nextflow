import smtplib
import os
import sys

from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email.utils import COMMASPACE, formatdate
from email import encoders

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--from-name')
parser.add_argument('-b', '--body', required=True)
parser.add_argument('--reply-to')
parser.add_argument('-t', '--to', action='append', nargs=1, metavar=('email'), required=True)
parser.add_argument('-s', '--subject', required=True)
parser.add_argument('-a', '--attach', action='append', nargs=1, metavar=('path'))
parser.add_argument('--body-start', default="")
args = parser.parse_args()

header = ["\n".join([x+',' for x in args.body_start.split(",")])]
footer = [
    "\n--\nThis message was sent automatically on behalf of the COVID-19 Genomics UK Consortium.",
    "If you believe you were not the intended recipient of this message, or if you no longer wish to receive these automated messages, please contact %s" % os.getenv("MAJE_MAINTAINER") + ", or reply to this message." if args.reply_to else ".",
]

# thanks https://stackoverflow.com/questions/3362600/how-to-send-email-attachments
msg = MIMEMultipart()
msg['To'] = ", ".join([x[0] for x in args.to])
msg['Date'] = formatdate(localtime=True)
msg['Subject'] = args.subject

if args.reply_to:
    msg["Reply-to"] = args.reply_to

if hasattr(args, "from_name"):
    msg["From"] = "(%s) <%s>" % (args.from_name, os.getenv("MAJE_USER"))
else:
    msg['From'] = os.getenv("MAJE_USER")

msg.attach(MIMEText("\n".join(header + [line.strip() for line in open(args.body).readlines()] + footer)))

for path in args.attach:
    path = path[0]
    part = MIMEBase('application', "octet-stream")
    print("attaching", path)
    with open(path, 'rb') as file:
        part.set_payload(file.read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition',
                    'attachment; filename="{}"'.format(os.path.basename(path)))
    msg.attach(part)

server = smtplib.SMTP_SSL('smtp.gmail.com')
server.ehlo()

server.login(os.getenv("MAJE_USER"), os.getenv("MAJE_PASS"))
server.sendmail(os.getenv("MAJE_USER"), [x[0] for x in args.to], msg.as_string())
server.close()
