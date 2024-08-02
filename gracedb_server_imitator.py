#!/usr/bin/env python

"""A small local HTTP server that pretends to be GraceDB and accepts event
uploads from PyCBC Live. At the moment this implement very minimal
functionality to let PyCBC Live upload without errors.

Start it on the same machine running PyCBC Live, and point PyCBC Live to
upload events to http://localhost:8000/api/."""

import logging
import http.server
import json
import re
from requests_toolbelt import MultipartDecoder


class FakeGraceDBServer(http.server.HTTPServer):
    """This class inherits from HTTPServer and keeps an internal list of G
    events that were uploaded to the server.
    """
    def __init__(self, name_port, handler):
        super().__init__(name_port, handler)
        self.events = []

    def new_event(self):
        gid = f"G{len(self.events)+1}"
        self.events.append(gid)
        logging.info("Created new event %s", gid)
        return gid


class MyHandler(http.server.BaseHTTPRequestHandler):
    def send_success(self):
        self.send_response(200)
        self.send_header('Content-type', 'application/json')
        self.end_headers()

    def do_GET(self):
        """Handle HTTP GET requests. Only /api/ is required.
        """
        output = None
        if self.path == "/api/":
            output = self.handle_api()
        if output is None:
            self.send_error(404)
        else:
            self.send_success()
            self.wfile.write(output.encode('utf-8'))

    def do_POST(self):
        """Handle HTTP POST requests. Only event uploads and log entries are
        required.
        """
        content_length = int(self.headers['Content-Length'])
        post_data = self.rfile.read(content_length)
        output = None
        if self.path == "/api/events/":
            # Create new event. The POST data is encoded in the complicated
            # "multipart form data" format, which includes the gzipped XML data
            # in the case of PyCBC Live uploads. This is just logged to stderr
            # for now, but one could decode the event data for added fun if
            # wanted.
            multipart_data = MultipartDecoder(post_data, self.headers['Content-Type'])
            for part in multipart_data.parts:
                logging.info('//// Multipart POST data for event upload ////')
                logging.info(part.headers)
                logging.info(part.content[:100])
            gid = self.server.new_event()
            # We need to return a JSON string with the GraceID for the client
            # to be happy.
            out_data = {
                "graceid": gid
            }
            output = json.dumps(out_data)
        elif re.search("/api/events/G.+/log", self.path):
            # Post a log entry (with an optional file upload).
            # Nothing happens at the moment.
            # One could save the content somewhere and actually check things.
            # We need to return a JSON string for things to work at the client.
            output = "{}"
        if output is None:
            self.send_error(404)
        else:
            self.send_success()
            self.wfile.write(output.encode('utf-8'))

    def handle_api(self):
        server_url_base = f"http://{self.server.server_name}:{self.server.server_port}/api/"
        api_versions = ["default", "v1", "v2"]
        data = {
            "links": {
                "events": server_url_base + "events/",
                "self": server_url_base,
                "performance": server_url_base + "performance/",
                "user-info": server_url_base + "user-info/"
            },
            "templates": {
                "event-detail-template": server_url_base + "events/{graceid}",
                "event-log-template": server_url_base + "events/{graceid}/log/",
                "event-log-detail-template": server_url_base + "events/{graceid}/log/{N}",
                "event-label-template": server_url_base + "events/{graceid}/labels/{label}",
            },
            "groups": [
                "CBC",
                "Test"
            ],
            "pipelines": [
                "pycbc"
            ],
            "searches": [
                "AllSky"
            ],
            "labels": [
                "EARLY_WARNING",
                "SNR_OPTIMIZED",
            ],
            "api-versions": api_versions,
            "server-version": "2.28.1",
            "API_VERSIONS": api_versions
        }
        return json.dumps(data)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S%z"
)

with FakeGraceDBServer(("", 8000), MyHandler) as httpd:
    httpd.serve_forever()
