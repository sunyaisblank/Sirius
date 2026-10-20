"""Read the real loader's returned float-control properties before dispatch."""
layout=json.loads(pathlib.Path(os.environ['SIRIUS_FLOAT_LAYOUT']).read_text())
property_rows=[]
property_fields=['denormBehaviorIndependence','roundingModeIndependence',
 'shaderSignedZeroInfNanPreserveFloat16','shaderSignedZeroInfNanPreserveFloat32','shaderSignedZeroInfNanPreserveFloat64',
 'shaderDenormPreserveFloat16','shaderDenormPreserveFloat32','shaderDenormPreserveFloat64',
 'shaderDenormFlushToZeroFloat16','shaderDenormFlushToZeroFloat32','shaderDenormFlushToZeroFloat64',
 'shaderRoundingModeRTEFloat16','shaderRoundingModeRTEFloat32','shaderRoundingModeRTEFloat64',
 'shaderRoundingModeRTZFloat16','shaderRoundingModeRTZFloat32','shaderRoundingModeRTZFloat64']
class PropertiesReturned(gdb.FinishBreakpoint):
    def __init__(self,pointer):
        super().__init__(gdb.newest_frame(),internal=True);self.pointer=pointer
    def stop(self):
        try:
            next_pointer=struct.unpack_from('<Q',memory(self.pointer,16),8)[0];seen=set();captured=[]
            while next_pointer:
                assert next_pointer not in seen and len(seen)<8;seen.add(next_pointer)
                header=memory(next_pointer,16);kind=struct.unpack_from('<I',header)[0]
                if kind==layout['sType']:
                    raw=memory(next_pointer,layout['bytes'])
                    values=dict(zip(property_fields,struct.unpack_from('<17I',raw,layout['first_field_offset'])))
                    assert values['shaderRoundingModeRTEFloat32']==values['shaderSignedZeroInfNanPreserveFloat32']==1
                    # This experiment never requests the unsupported denormal mode.
                    assert values['shaderDenormPreserveFloat32']==0
                    captured.append({'pointer':next_pointer,'raw_hex':raw.hex(),'values':values})
                next_pointer=struct.unpack_from('<Q',header,8)[0]
            assert len(captured)==1
            property_rows.extend(captured)
            (feature_folder/'float-controls-properties.json').write_text(json.dumps(property_rows,indent=2)+'\n')
            return False
        except Exception as error:return failure(error)
class PropertiesQuery(gdb.Breakpoint):
    def stop(self):
        try:
            feature_boundary();PropertiesReturned(reg('rsi'));return False
        except Exception as error:return failure(error)
PropertiesQuery('vkGetPhysicalDeviceProperties2',internal=True)
