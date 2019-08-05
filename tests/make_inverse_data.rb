require 'matrix'

def print_c_source1(m)
  print "    {\n"
  print "      #{m.row_size},\n"
  print "      (double[]) {\n"

  m.to_a.each {|row|
    tmp = row.inject([]) {|m, n| m << ("% 4d" % n)}
    print "        #{tmp.join(",")},\n"
  }

  print "      },\n"
  print "    },\n"
end

def print_c_source2(m)
  print "    {\n"
  print "      #{m.row_size},\n"
  print "      (double[]) {\n"

  m.to_a.each {|row|
    tmp = row.inject([]) {|m, n| m << ("% 34.30f" % n)}
    print "        #{tmp.join(",")},\n"
  }

  print "      },\n"
  print "    },\n"
end

print <<~EOT
  typedef struct {
    int size;
    double* val;
  } matrix_info_t;

  static struct {
    matrix_info_t op;
    matrix_info_t ans;
  } data[] = {
EOT

2000.times {

  begin
    n  = rand(15) + 1
    op = Matrix[*(Array.new(n) {Array.new(n) {rand(-10...+10)}})]
  end until op.regular?

  ans = op.inverse

  print "  {\n"
  print_c_source1(op)
  print_c_source2(ans)
  print "  },\n"
}

print <<~EOT
 };
EOT
