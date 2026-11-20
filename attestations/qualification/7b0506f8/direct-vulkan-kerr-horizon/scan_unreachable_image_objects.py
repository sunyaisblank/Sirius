"""Inspect local Git blobs outside all reachable refs for original-tape leads."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys

work = Path(sys.argv[1])
all_objects = work/'all-objects.tsv'
reachable_objects = work/'reachable.txt'
reachable = set(reachable_objects.read_text().splitlines())
inventory = [tuple(line.split()) for line in all_objects.read_text().splitlines()]
unreachable = sorted((identity,int(size)) for kind,identity,size in inventory
                     if kind=='blob' and identity not in reachable)
assert len(unreachable) <= 100 and sum(sz for _,sz in unreachable) <= 20*1024*1024
matches = []
counts = Counter()
signatures = {'png':b'\x89PNG\r\n\x1a\n', 'openexr':bytes.fromhex('762f3101'),
              'gzip':b'\x1f\x8b', 'zip':b'PK\x03\x04',
              '7zip':bytes.fromhex('377abcaf271c')}
terms = (b'rebuild-baseline', b'kerr-png.png', b'kerr-exr.exr', b'alcubierre-png.png',
         b'original tape', b'pre-rebuild')
for identity,size in unreachable:
    raw = subprocess.check_output(['git','cat-file','blob',identity])
    assert len(raw)==size
    storage = [key for key,marker in signatures.items() if raw.startswith(marker)]
    embedded = [key for key,marker in signatures.items()
                if key in {'png','openexr'} and raw.find(marker,1)!=-1]
    text_markers = []
    try:
        decoded = raw.decode('utf-8')
        counts['utf8']+=1
        lower = decoded.lower().encode()
        text_markers = [term.decode() for term in terms if term in lower]
    except UnicodeDecodeError:
        counts['binary']+=1
    if storage or embedded or text_markers:
        matches.append({'blob':identity,'bytes':size,'format':storage,'embedded_signatures':embedded,
                        'forensic_text_terms':text_markers})
report = {'issue':100,
          'revision':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
          'scope':'All available Git blob objects outside refs reachable from git rev-list --all. '
                  'Magic identifies a direct container lead. Text terms and embedded matches require inspection.',
          'all_object_inventory_sha256':hashlib.sha256(all_objects.read_bytes()).hexdigest(),
          'reachable_objects':len(reachable),
          'unreachable_blob_count':len(unreachable),'uncompressed_blob_bytes':sum(sz for _,sz in unreachable),
          'unreachable_inventory_sha256':hashlib.sha256(''.join(f'{o} {n}\n' for o,n in unreachable).encode()).hexdigest(),
          'types':counts,'candidates':matches}
(work/'scan.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print(json.dumps(report,sort_keys=True,indent=2))
