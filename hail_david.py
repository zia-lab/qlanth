#!/usr/bin/env python3

import http.client, urllib
import argparse

try:
  from dave import secrets
  push_token = secrets['pushover_token']
  push_user  = secrets['pushover_user']
except:
  print("Ask David about secrets.")

def send_message(message, app='Ringer'):
    '''
    Send a notification through Pushover.

    Parameters
    ----------
    message : str
        The message to be sent.
    app : str
        The app to which the message will be sent. Default is Ringer.
    
    Returns
    None
    '''
    token_dict = {'Ringer': push_token}
    app_token = token_dict[app]
    conn = http.client.HTTPSConnection("api.pushover.net",443)
    endpoint = "/1/messages.json"
    conn.request("POST", endpoint,
      urllib.parse.urlencode({
        "token": app_token,
        "user": push_user,
        "message": message,
      }), { "Content-type": "application/x-www-form-urlencoded" })
    return conn.getresponse().read().decode()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Send me a message.')
    parser.add_argument('message', type=str, help='Message to be sent.')
    args = parser.parse_args()
    send_message(args.message)