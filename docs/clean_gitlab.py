#!/usr/bin/env python3

''' Clean pipelines in Gitlab
'''

import time

import requests

project_id = '12606413'
token = 'glpat-md4-1yPptUs33w7PogPY'
server = 'gitlab.com'

print("Creating list of all jobs that currently have artifacts...")
# We skip the first page.
url = f"https://{server}/api/v4/projects/{project_id}/pipelines?per_page=100&page=1"
while url:
    print(f"Processing page: {url}")
    response = requests.get(
        url,
        headers={
            'private-token': token,
        },
    )

    if response.status_code in [500, 429]:
        print(f"Status {response.status_code}, retrying.")
        time.sleep(10)
        continue

    response.raise_for_status()
    response_json = response.json()
    for job in response_json[1:]:
        print(job)
        #if job.get('id', None) != 408045050:
        job_id = job['id']
        delete_response = requests.delete(
            f"https://{server}/api/v4/projects/{project_id}/pipelines/{job_id}",
            headers={
                'private-token': token,
            },
        )
        print(f"Processing job ID: {job_id} - status: {delete_response.status_code}")

    url = response.links.get('next', {}).get('url', None)
