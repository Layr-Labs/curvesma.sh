/* layout_assert.c -- the layout law, checked at COMPILE time.
 * Included by scripts/check_layout.py: if rec_t ever diverges from the stride
 * and offsets the listing evidence/step_secp_v15g.s walks, this TU fails to
 * compile, so the package cannot be regenerated. */
#include <stdio.h>
#include "rec_layout.h"
int main(void) {
    printf("{\"sizeof_rec_t\": %zu, \"off_x\": %zu, \"off_y\": %zu, \"off_idx\": %zu, "
           "\"off_len\": %zu, \"off_esc\": %zu}\n",
           sizeof(rec_t), __builtin_offsetof(rec_t, x), __builtin_offsetof(rec_t, y),
           __builtin_offsetof(rec_t, idx), __builtin_offsetof(rec_t, len),
           __builtin_offsetof(rec_t, esc));
    return 0;
}
