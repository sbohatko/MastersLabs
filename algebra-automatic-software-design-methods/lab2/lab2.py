class Mealy(object):
    """Mealy Machine : Finite Automata with Output """

    def __init__(self, states, input_alphabet, output_alphabet, transitions, initial_state):
        """
        6 tuple (Q, ∑, O, δ, X, q0) where −

        states is a finite set of states.

        alphabet is a finite set of symbols called the input alphabet.

        output_alphabet is a finite set of symbols called the output alphabet.

        transitions is the resultant data dictionary of input and output transition functions

        initial_state is the initial state from where any input is processed (q0 ∈ Q).
        """
        self.states = states
        self.input_alphabet = input_alphabet
        self.output_alphabet = output_alphabet
        self.transitions = transitions
        self.initial_state = initial_state

    def get_output_from_string(self, string):
        """Return Mealy Machine's output when a given string is given as input"""

        temp_list = list(string)
        current_state = self.initial_state
        output = ''
        for x in temp_list:
            output+= self.transitions[current_state][x][1]
            current_state = self.transitions[current_state][x][0]

        return output

    def convert_to_moore(self):
        moore_transitions = {}
        temp_list = []
        moore_output_table = {}
        moore_initial_state = self.initial_state
        for x in self.transitions.keys():
            for a in self.input_alphabet:
                temp_list.append(self.transitions[x][a])

        temp_list_2 = []
        for x in temp_list:
            for y in temp_list:
                if x[0] == y[0] and x[1] != y[1]:
                    if x not in temp_list_2 and y not in temp_list_2:
                        temp_list_2.append(x)
                        temp_list_2.append(y)

        temp_list_3 = []
        for x in temp_list_2:
            if x[0] not in temp_list_3:
                temp_list_3.append(x[0])

        if self.initial_state in temp_list_3:
            moore_initial_state = self.initial_state + self.output_alphabet[0]

        for x in temp_list_2:
            for a in self.input_alphabet:
                if self.transitions[x[0]][a][0] in temp_list_3:
                    next_state = self.transitions[x[0]][a][0]
                    output = self.transitions[x[0]][a][1]

                    next_state = next_state + output
                    try:
                        moore_transitions[x[0] + x[1]][a] = next_state
                    except KeyError as e:
                        moore_transitions[x[0] + x[1]] = {}
                        moore_transitions[x[0] + x[1]][a] = next_state

                    if next_state not in moore_output_table.keys():
                        moore_output_table[next_state] = output

                else:
                    try:
                        moore_transitions[x[0] + x[1]][a] = self.transitions[x[0]][a][0]
                    except KeyError as e:
                        moore_transitions[x[0] + x[1]] = {}
                        moore_transitions[x[0] + x[1]][a] = self.transitions[x[0]][a][0]

                    if moore_transitions[x[0] + x[1]][a] not in moore_output_table.keys():
                        moore_output_table[moore_transitions[x[0] + x[1]][a]] = self.transitions[x[0]][a][1]

        for x in self.transitions.keys():
            if x not in moore_transitions.keys() and x not in temp_list_3:
                for a in self.input_alphabet:
                    if self.transitions[x][a][0] in temp_list_3:
                        next_state = self.transitions[x][a][0]
                        output = self.transitions[x][a][1]

                        next_state = next_state + output
                        try:
                            moore_transitions[x][a] = next_state
                        except KeyError as e:
                            moore_transitions[x] = {}
                            moore_transitions[x][a] = next_state

                        if next_state not in moore_output_table.keys():
                            moore_output_table[next_state] = output

                    else:
                        try:
                            moore_transitions[x][a] = self.transitions[x][a][0]
                        except KeyError as e:
                            moore_transitions[x] = {}
                            moore_transitions[x][a] = self.transitions[x][a][0]

                        if self.transitions[x][a][0] not in moore_output_table.keys():
                            moore_output_table[self.transitions[x][a][0]] = self.transitions[x][a][1]

        moore_states = []
        for s in moore_transitions.keys():
            if s not in moore_states:
                moore_states.append(s)



        moore_from_mealy = Moore(
            moore_states,
            self.input_alphabet,
            self.output_alphabet,
            moore_transitions,
            moore_output_table,
            moore_initial_state
        )

        print(moore_from_mealy)


    def __str__(self):
        output = "\nMealy Machine" + \
                 "\nStates " + str(self.states) + \
                 "\nTransitions " + str(self.transitions) + \
                 "\nInital State " + str(self.initial_state) + \
                 "\nInital Alphabet " + str(self.input_alphabet) + \
                 "\nOutput Alphabet" + str(self.output_alphabet)

        return output

#--------------------------------------------------------------------------------------

class Moore(object):
    """Moore Machine : Finite Automata with Output"""

    def __init__(self, states, input_alphabet, output_alphabet, transitions, initial_state, output_table ):
        """
        states: Finite set of states
        input_alphabet: Alphabet of letters for forming input string
        output_alphabet: Alphabet of letters for forming output characters
        transitions: Transition Table
        output_table: Output Table to show what character from output_alphabet is printed when a state from 'states'
        is reached
        """

        self.states = states
        self.input_alphabet = input_alphabet
        self.output_alphabet = output_alphabet
        self.transitions = transitions
        self.output_table = output_table
        self.initial_state = initial_state


    def __str__(self):
        """"Pretty Print the Moore Machine"""

        output = "\nMoore Machine" + \
                 "\nStates = " + str(self.states) + \
                 "\nInput Alphabet = " + str(self.input_alphabet) + \
                 "\nOutput Alphabet = " + str(self.output_alphabet) + \
                 "\nTransitions = " + str(self.transitions) + \
                 "\nInitial State = " + str(self.initial_state) + \
                 "\nOutput Table = " + str(self.output_table)

        return output


mealy_2 = Mealy(['a','b'],
                ['x', 'y'],
                ['0', '1'],
                {
                    'a' : {
                        'x' : ('b', '0'),
                        'y' : ('a', '1')
                    },
                    'b' : {
                        'x' : ('a', '0'),
                        'y' : ('b', '1')
                    }
                },
                'a'
                )
print(mealy_2)
print(mealy_2.convert_to_moore())