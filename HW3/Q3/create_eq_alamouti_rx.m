function H_eq = create_eq_alamouti_rx(H)
H_eq = [H(1),H(2);
        H(2)',-H(1)'];
end