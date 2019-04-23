#!/usr/bin/env python

"""Walk through a directory of PyCBC Live coinc triggers and make their
coinc_event_ids unique."""

import sys
import os
import glob
import gzip
import tqdm


in_fn = glob.glob(sys.argv[1] + '/*.xml*')

for i, fn in enumerate(tqdm.tqdm(in_fn)):
    with gzip.open(fn, 'rb') as f_in:
        text = f_in.read()

    if 'coinc_event:coinc_event_id:0' not in text:
        continue

    fixed_text = text.replace('coinc_event:coinc_event_id:0',
                              'coinc_event:coinc_event_id:{}'.format(i))

    fn_out = os.path.join(sys.argv[2], os.path.basename(fn))
    with gzip.open(fn_out, 'wb') as f_out:
        f_out.write(fixed_text)
