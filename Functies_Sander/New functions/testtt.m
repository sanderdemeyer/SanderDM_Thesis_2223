test1(4, 'test', 5)
test2(5, 'D1', 7, 'test', 11)

function zz = test1(a, kwargs)
    arguments
        a
        kwargs.D1 = 0
        kwargs.D2 = 0
        kwargs.test = 9
    end
    zz = a + kwargs.test;
end

function zzz = test2(b, kwargs)
    arguments
        b
        kwargs.D1 = 0
        kwargs.test = 9
    end
    disp(kwargs);
    disp(kwargs{:});
    zzz = test1(b+1, kwargs{:});
end