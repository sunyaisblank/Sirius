"""Read the original CPU writer test's live files in place before its cleanup."""

import ctypes
import hashlib
import importlib.util
import json
from math import isfinite
import os
from pathlib import Path
import select
import struct
import sys
import time

ROOT=Path.cwd()
WORK=ROOT/'out/strict-backend-science'
OWNER=WORK/'ordinary/host-controls/runs/host-controls/owner.json'
RENDERS=ROOT/'renders/strict-backend-science/ordinary/host-controls/host-controls'
PNG_AUTHORITY=ROOT/'scripts/verify-attestation.py'
EXR_AUTHORITY=ROOT/'out/science-replay/audit_exr.py'


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pin(raw):
    return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}


def decode_exr(raw,helper):
    assert struct.unpack_from('<II',raw)==(20000630,2)
    offset=8
    fields={}
    while raw[offset]:
        name,offset=helper.field_name(raw,offset)
        kind,offset=helper.field_name(raw,offset)
        size,=struct.unpack_from('<I',raw,offset)
        offset+=4
        assert name not in fields and offset+size<=len(raw)
        fields[name]=(kind,raw[offset:offset+size])
        offset+=size
    offset+=1
    assert fields['dataWindow']==('box2i',struct.pack('<4i',0,0,63,63))
    assert fields['compression']==('compression',b'\x03')
    assert fields['lineOrder'][1]==b'\x00'
    kind,chlist=fields['channels']
    assert kind=='chlist'
    names=[]
    at=0
    while chlist[at]:
        name,at=helper.field_name(chlist,at)
        assert struct.unpack_from('<i4xii',chlist,at)==(2,1,1)
        names.append(name)
        at+=16
    assert sorted(names)==list('ABGR') and at+1==len(chlist)
    pixels={name:[] for name in names}
    blocks=struct.unpack_from('<4Q',raw,offset)
    next_block=offset+32
    for block,position in enumerate(blocks):
        assert position==next_block
        row,packed=struct.unpack_from('<ii',raw,position)
        assert row==block*16 and packed>0
        next_block=position+8+packed
        values=helper.uncompress(raw[position+8:next_block],16*64*4*4)
        floats=struct.unpack('<4096f',values)
        for index in range(16):
            for channel,name in enumerate(names):
                at=(index*4+channel)*64
                pixels[name].extend(floats[at:at+64])
    assert next_block==len(raw)
    assert all(len(values)==4096 and all(map(isfinite,values)) for values in pixels.values())
    rgb=[value for name in 'RGB' for value in pixels[name]]
    assert max(rgb)>min(rgb) and max(rgb)>0 and all(value==1 for value in pixels['A'])
    return dict(width=64,height=64,channels=names,format='ZIP scanline RGBA IEEE binary32',
                finite_samples=16384,rgb_min=min(rgb),rgb_max=max(rgb),
                distinct_rgb_triplets=len(set(zip(*(pixels[name] for name in 'RGB')))))


def main():
    png=load(PNG_AUTHORITY,'qualified_png_reader')
    exr=load(EXR_AUTHORITY,'qualified_exr_reader')
    result={'sources':{
        str(p.relative_to(ROOT)):pin(p.read_bytes()) for p in (PNG_AUTHORITY,EXR_AUTHORITY)},
        'producer':'RenderSessionProbe.CpuKerrRenderProducesValidPngAndExr','captures':{}}
    destination=WORK/'writer-audit.json'

    def save(kind,path,raw,data):
        assert len(raw)<200_000
        result['captures'][kind]={'path':str(path.relative_to(ROOT)), 'file':pin(raw),'decode':data}
        destination.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
        print(f'{kind}: {pin(raw)}',flush=True)

    image,=RENDERS.glob('sirius-render-probe-*/probe_kerr.png')
    raw=image.read_bytes()
    png_data=png.inspect_png(image)
    assert png_data['dimensions']==(64,64) and png_data['distinct_colours']>1
    assert raw==image.read_bytes()
    save('png',image,raw,png_data)

    # Retain a read-only descriptor on creation. The test can delete its
    # temporary filename immediately after its own decoder has finished.
    create,closed,moved=0x100,0x08,0x80
    target=image.with_suffix('.exr')
    linux=ctypes.CDLL(None,use_errno=True)
    notification=linux.inotify_init1(os.O_CLOEXEC)
    if notification<0:
        raise OSError(ctypes.get_errno(),'Cannot watch the original writer')
    reader=None
    try:
        watch=linux.inotify_add_watch(notification,os.fsencode(target.parent),create|closed|moved)
        if watch<0:
            raise OSError(ctypes.get_errno(),'Cannot watch the original output directory')
        start=time.monotonic()
        while time.monotonic()-start<1400:
            ready,_,_=select.select([notification],[],[],.1)
            if ready:
                data=os.read(notification,8192)
                offset=0
                while offset<len(data):
                    _,mask,_,length=struct.unpack_from('iIII',data,offset)
                    offset+=16
                    name=data[offset:offset+length].rstrip(b'\0')
                    offset+=length
                    if name!=b'probe_kerr.exr':
                        continue
                    if reader is None:
                        try:
                            reader=target.open('rb')
                        except FileNotFoundError:
                            continue
                    if mask & (closed|moved):
                        payload=reader.read(200_001)
                        save('exr',target,payload,decode_exr(payload,exr))
                        return 0
            if json.loads(OWNER.read_text())['status']=='terminal':
                return 1
        return 1
    finally:
        if reader is not None:
            reader.close()
        os.close(notification)


if __name__=='__main__':
    sys.exit(main())
