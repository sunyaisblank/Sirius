"""Bind this fresh cohort to the installed zlib update; preserve old receipts."""
from pathlib import Path
import hashlib
import json
import re
import subprocess


def apply(root, work, expected):
    name = '/usr/lib/x86_64-linux-gnu/libz.so.1.3'
    old = {'bytes': 113000, 'sha256': '86200da370f20476a2507e9097a789b5ef97269b4ca8d5e164ad82dab9d99892'}
    current = {'bytes': 113000, 'sha256': 'f61677c7bdb7c2aa2a243ea781bfd16afcb361cede7b4ef76127a9d8eb9a5242'}
    def pin(path):
        with path.open('rb') as stream:
            return {'bytes': path.stat().st_size, 'sha256': hashlib.file_digest(stream, 'sha256').hexdigest()}
    incoming = expected[name]
    assert incoming in (old, current) and pin(Path(name)) == current
    cpp_preparation = root / 'out/native-cpu-arithmetic-closure/prepared-cpp.json'
    cpp = json.loads(cpp_preparation.read_text())
    assert cpp['inputs'][name] == current
    version = subprocess.check_output(['dpkg-query', '-W', '-f=${Version}', 'zlib1g:amd64'], text=True)
    assert version == '1:1.3.dfsg-3.1ubuntu2.3'
    verification = subprocess.check_output(['dpkg', '-V', 'zlib1g:amd64'], text=True)
    assert not verification.strip()
    log = [line for line in Path('/var/log/dpkg.log').read_text().splitlines()
           if '2026-10-25 06:06:09 upgrade zlib1g:amd64 ' in line]
    assert log == ['2026-10-25 06:06:09 upgrade zlib1g:amd64 1:1.3.dfsg-3.1ubuntu2.2 1:1.3.dfsg-3.1ubuntu2.3']
    host = json.loads((root / 'out/native-cpu-upload-aware-staging/host-build-result.json').read_text())
    provider = json.loads((root / 'out/native-cpu-arithmetic-closure/provider/provider-build-result.json').read_text())
    baseline = json.loads((root / 'out/native-cpu-arithmetic-sharing/providers.json').read_text())['providers']['baseline']
    products = [(root / host['executable'], host['executable_identity']),
                (root / host['frame_executable'], host['frame_executable_identity']),
                (root / provider['library'], provider['library_identity']),
                (Path(baseline['library']), baseline['library_identity'])]
    resolutions = []
    for path, identity in products:
        assert pin(path) == identity
        output = subprocess.check_output(['/usr/bin/ldd', '-r', str(path)], text=True, stderr=subprocess.STDOUT)
        assert 'undefined symbol' not in output and 'not found' not in output
        resolutions.append({'path': str(path), 'identity': identity,
                            'resolution': re.sub(r'0x[0-9a-fA-F]+', '0xADDRESS', output)})
    assert pin(Path(name)) == current
    record = {'path': name, 'historical_build_dependency': old, 'current_runtime': current,
              'incoming_merged_dependency': incoming, 'current_CPP_build_preparation': pin(cpp_preparation),
              'package': 'zlib1g:amd64', 'version': version, 'package_integrity_output': verification,
              'upgrade_log': log, 'existing_product_resolution': resolutions,
              'scope': 'Historical P host receipts pin the old library; the new dense C++ build explicitly pins the current installed library. Preserve both provenances and verify the incoming merged pin against those exact two identities. This cohort uses the current package after product/library-resolution checks; no package change is performed. All other identities, original C8/full38/frame assertions/counts/governors remain. Z broader gates are predecessor supporting evidence only. No historical qualification, numerical equivalence or release claim is transferred across the package update.'}
    with (work / 'runtime-transition.json').open('x') as stream:
        json.dump(record, stream, indent=2)
        stream.write('\n')
    expected[name] = current
