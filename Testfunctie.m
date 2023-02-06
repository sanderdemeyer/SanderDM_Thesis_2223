start_matrix = [1 1 1 1; 0 0 0 0; 0 0 0 0; 0 0 0 0];
disp(test_fun(start_matrix))


fun = @root2d;
disp(fun(5));

a = fsolve(test_fun, start_matrix);

function a = test_fun(x)
    disp(x);
    a = x - transpose(x);
end