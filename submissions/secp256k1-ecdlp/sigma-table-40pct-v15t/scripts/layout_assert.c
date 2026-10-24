/* layout_assert.c -- the layout law, checked at COMPILE time.
 * Included by scripts/check_layout.py: if rec_t ever diverges from the stride
 * and offsets the listing evidence/step_secp_v15g.s walks, this TU fails to
 * compile, so the package cannot be regenerated. */
#include <stdio.h>
#include "rec_layout.h"
int main(void) {
    printf("{\"sizeof_rec_t\": %zu, \"off_x\": %zu, \"off_y\": %zu, \"off_idx\": %zu, "
           "\"off_len\": %zu, \"off_esc\": %zu, "
           "\"sizeof_worker_t\": %zu, \"worker_off_rec\": %zu, \"worker_off_a\": %zu, "
           "\"worker_off_c\": %zu, \"worker_off_prev0\": %zu, \"worker_off_prev1\": %zu, "
           "\"worker_off_minx\": %zu, \"worker_off_rep\": %zu, \"worker_off_win\": %zu, "
           "\"worker_batch\": %u, \"prefix_bytes\": %u}\n",
           sizeof(rec_t), __builtin_offsetof(rec_t, x), __builtin_offsetof(rec_t, y),
           __builtin_offsetof(rec_t, idx), __builtin_offsetof(rec_t, len),
           __builtin_offsetof(rec_t, esc),
           sizeof(worker_t), __builtin_offsetof(worker_t, rec),
           __builtin_offsetof(worker_t, a), __builtin_offsetof(worker_t, c),
           __builtin_offsetof(worker_t, prev0), __builtin_offsetof(worker_t, prev1),
           __builtin_offsetof(worker_t, minx), __builtin_offsetof(worker_t, rep),
           __builtin_offsetof(worker_t, win), WORKER_BATCH, PREFIX_BYTES);
    return 0;
}
