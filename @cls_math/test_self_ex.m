function result = test_self_ex(this, num_in)
    result = num_in + 999;
    fprintf('%s says %d\n', this.name, result);
end

