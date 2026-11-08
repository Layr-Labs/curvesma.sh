import json
import sys
import random

def floyd_rho(order, seed=42):
    random.seed(seed)
    # Additive cyclic group; recover discrete log of a random target
    target = random.randrange(1, order)
    def step(x, c, d):
        # Simple illustrative step (not a real partition)
        return (x + 1) % order, (c + 1) % order, d
    # Tortoise and hare
    t_x, t_c, t_d = 0, 0, 0
    h_x, h_c, h_d = 0, 0, 0
    steps = 0
    max_steps = int(8 * (order ** 0.5)) + 20
    while steps < max_steps:
        # tortoise one step
        t_x, t_c, t_d = step(t_x, t_c, t_d)
        # hare two steps
        h_x, h_c, h_d = step(h_x, h_c, h_d)
        h_x, h_c, h_d = step(h_x, h_c, h_d)
        steps += 1
        if t_x == h_x and steps > 1:
            # Cycle detected; in a real attack we would solve for the log
            # Here we just confirm detection
            return {"success": True, "steps": steps, "order": order, "target": target}
    return {"success": False, "steps": steps, "order": order, "target": target}

if __name__ == "__main__":
    data = json.load(sys.stdin)
    order = data["order"]
    result = floyd_rho(order)
    print(json.dumps(result))
