"""Observe controls, plus the first real compute callback's state restoration."""
layout=json.loads(pathlib.Path(os.environ['SIRIUS_FLOAT_LAYOUT']).read_text())
expected_denorm=int(os.environ['SIRIUS_EXPECT_DENORM'])
property_fields=['denormBehaviorIndependence','roundingModeIndependence',
 'shaderSignedZeroInfNanPreserveFloat16','shaderSignedZeroInfNanPreserveFloat32','shaderSignedZeroInfNanPreserveFloat64',
 'shaderDenormPreserveFloat16','shaderDenormPreserveFloat32','shaderDenormPreserveFloat64',
 'shaderDenormFlushToZeroFloat16','shaderDenormFlushToZeroFloat32','shaderDenormFlushToZeroFloat64',
 'shaderRoundingModeRTEFloat16','shaderRoundingModeRTEFloat32','shaderRoundingModeRTEFloat64',
 'shaderRoundingModeRTZFloat16','shaderRoundingModeRTZFloat32','shaderRoundingModeRTZFloat64']
properties=[];callbacks=[];seen_workers=set()
def save_states():
    (feature_folder/'compute-boundaries.json').write_text(json.dumps(callbacks,indent=2)+'\n')
class PropertiesReturned(gdb.FinishBreakpoint):
    def __init__(self,pointer):super().__init__(gdb.newest_frame(),internal=True);self.pointer=pointer
    def stop(self):
        try:
            next_pointer=struct.unpack_from('<Q',memory(self.pointer,16),8)[0];seen=set();found=[]
            while next_pointer:
                assert next_pointer not in seen and len(seen)<8;seen.add(next_pointer);header=memory(next_pointer,16)
                if struct.unpack_from('<I',header)[0]==layout['sType']:
                    raw=memory(next_pointer,layout['bytes']);values=dict(zip(property_fields,struct.unpack_from('<17I',raw,layout['first_field_offset'])))
                    assert values['shaderDenormPreserveFloat32']==expected_denorm
                    assert values['shaderRoundingModeRTEFloat32']==values['shaderSignedZeroInfNanPreserveFloat32']==1
                    found.append({'raw_hex':raw.hex(),'values':values})
                next_pointer=struct.unpack_from('<Q',header,8)[0]
            assert len(found)==1;properties.extend(found);(feature_folder/'float-controls-properties.json').write_text(json.dumps(properties,indent=2)+'\n');return False
        except Exception as error:return failure(error)
class PropertiesQuery(gdb.Breakpoint):
    def stop(self):
        try:feature_boundary();PropertiesReturned(reg('rsi'));return False
        except Exception as error:return failure(error)
class ComputeReturned(gdb.FinishBreakpoint):
    def __init__(self,row):super().__init__(gdb.newest_frame(),internal=True);self.row=row
    def stop(self):
        try:
            self.row.update(returned=True,mxcsr_after=reg('mxcsr'));save_states();return False
        except Exception as error:return failure(error)
class ComputeEntered(gdb.Breakpoint):
    def stop(self):
        try:
            worker=gdb.selected_thread().global_num
            if worker not in seen_workers:
                seen_workers.add(worker);library=gdb.solib_name(reg('rip'))
                assert library and pathlib.Path(library).resolve(strict=True)==pathlib.Path(os.environ['SIRIUS_ACTUAL_PROVIDER']).resolve(strict=True)
                row={'worker':worker,'mxcsr_before':reg('mxcsr'),'returned':False};callbacks.append(row);ComputeReturned(row);save_states()
            return False
        except Exception as error:return failure(error)
PropertiesQuery('vkGetPhysicalDeviceProperties2',internal=True)
ComputeEntered('cs_exec_fn',internal=True)
