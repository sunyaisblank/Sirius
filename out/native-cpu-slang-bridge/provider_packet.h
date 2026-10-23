/* Included by lp_state_cs.c after lp_cs_job_info is defined. */
#include "sirius_native_bridge.h"

static bool
lp_native_cpu_packet(const struct lp_cs_job_info *job,
                     uint32_t x, uint32_t y, uint32_t z)
{
   if (!job->native_cpu_stage || job->native_cpu_stage > 3 ||
       job->block_size[0] != 64 || job->block_size[1] != 1 || job->block_size[2] != 1 ||
       job->req_local_mem || job->io || job->payload ||
       y || z || job->grid_size[1] != 1 || job->grid_size[2] != 1 ||
       x >= job->grid_size[0])
      return false;

   /* Vulkan set0 is constants[1], not the Gallium SSBO array. This screen's
    * lp_get_constant_buffer_stride is exactly sizeof(float). The admitted
    * layout has exactly two descriptors and no inline uniform block.
    */
   const struct lp_jit_buffer *table = &job->current->jit_resources.constants[1];
   _Static_assert(2 * sizeof(struct lp_descriptor) % sizeof(float) == 0,
                  "whole descriptor-table stride");
   const size_t table_bytes = 2 * sizeof(struct lp_descriptor);
   if (!table->u || table->num_elements != table_bytes / sizeof(float))
      return false;
   uintptr_t table_address = (uintptr_t)table->u;
   if (table_address > UINTPTR_MAX - table_bytes)
      return false;

   struct lp_descriptor descriptors[2];
   memcpy(descriptors, table->u, sizeof(descriptors));
   const struct lp_jit_buffer *input = &descriptors[0].buffer;
   const struct lp_jit_buffer *output = &descriptors[1].buffer;
   uintptr_t output_address = (uintptr_t)output->u;
   if (output_address > UINTPTR_MAX - output->num_elements ||
       (output_address < table_address + table_bytes &&
        table_address < output_address + output->num_elements))
      return false;
   /* SSBO descriptor extents are bytes. Full numerical/program/range checks
    * occur together in try_packet before any output write. A refusal never
    * means that an earlier row was executed; the ordinary JIT remains safe.
    */
   uint32_t groups = job->intergroup_row_owner ? MIN2(8, job->grid_size[0] - x) : 1;
   switch (job->native_cpu_stage) {
   case 1:
      return sirius_native_transport_try_packet(input->u, input->num_elements,
         (void *)output->u, output->num_elements, x, y, z,
         job->grid_size[0], job->grid_size[1], job->grid_size[2], groups);
   case 2:
      return sirius_native_endpoint_try_packet(input->u, input->num_elements,
         (void *)output->u, output->num_elements, x, y, z,
         job->grid_size[0], job->grid_size[1], job->grid_size[2], groups);
   case 3:
      return sirius_native_dopri_phase_try_packet(input->u, input->num_elements,
         (void *)output->u, output->num_elements, x, y, z,
         job->grid_size[0], job->grid_size[1], job->grid_size[2], groups);
   }
   return false;
}
