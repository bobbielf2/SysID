function y = myModel(th,test_case)

switch test_case
    case 1
        y = sin(th);
    case 2
        U = reactDiffuse1d([th(1), th(1), th(2)]);
        y = U(:);
    case 3
        U = reactDiffuse1d(th);
        y = U(:);
end