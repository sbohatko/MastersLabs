def epsilon_closure(states, transitions, epsilon='e'):
    e_closure = set(states)

    stack = []
    for state in states:
        stack.append(state)

    while stack:
        state = stack.pop()
        key = (state, epsilon)
        if key in transitions:
            for s in transitions[key]:
                if s not in e_closure:
                    e_closure.add(s)
                    stack.append(s)

    return tuple(sorted(list(e_closure)))


def nfa_to_dfa(nfa):
    nfa_states = nfa["A"]
    alphabet = nfa["X"]
    transitions = nfa["f"]
    start_state = nfa["a0"]
    nfa_final_states = nfa["F"]

    dfa_states = []
    dfa_transitions = {}
    unprocessed_states = []

    start_state_e_closure = epsilon_closure([start_state], transitions)

    dfa_states.append(start_state_e_closure)
    unprocessed_states.append(start_state_e_closure)

    while unprocessed_states:
        current_state = unprocessed_states.pop(0)

        for symbol in alphabet:
            next_states = set()

            for state in current_state:
                if (state, symbol) in transitions:
                    for s in transitions[(state, symbol)]:
                        next_states.add(s)

            next_states_e_closure = epsilon_closure(next_states, transitions)

            if next_states_e_closure not in dfa_states and len(next_states_e_closure) > 0:
                dfa_states.append(next_states_e_closure)
                unprocessed_states.append(next_states_e_closure)

            dfa_transitions[(current_state, symbol)] = next_states_e_closure

    dfa_final_states = [state for state in dfa_states if any(s in nfa_final_states for s in state)]

    return {
        "A": dfa_states,
        "X": alphabet,
        "f": dfa_transitions,
        "a0": start_state_e_closure,
        "F": dfa_final_states
    }


# Given NFA
nfa = {"A": ["1", "2", "3", "4"], "X": ["a", "b", "c"], "f": {('1', 'b'): ["2", "3"], ('2', 'a'): ["2"], ('3', 'c'): ["4"], ('4', 'a'): ["4"], ('1', 'a'): ["2"], ('1', 'c'): ["4"], ('2', 'b'): ["2", "3"]}, "a0": "1", "F": ["1", "3", "4"], "epsilon": "e"}


# Convert NFA to DFA
dfa = nfa_to_dfa(nfa)
print("NFA: ", nfa)
print("\nDFA: ", dfa)

# Print DFA
print("\nDFA States:")
for state in dfa["A"]:
    print(state)

print("\nAlphabet (X):")
for symbol in dfa["X"]:
    print(symbol)

print("\nTransitions (f):")
for transition, next_state in dfa["f"].items():
    print(f"Transition from {transition[0]} with input {transition[1]} goes to {next_state}")

print("\nStart State (a0):")
print(dfa["a0"])

print("\nFinal States (F):")
for state in dfa["F"]:
    print(state)