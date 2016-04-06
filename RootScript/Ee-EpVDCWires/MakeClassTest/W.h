
class W{
public:
  W(){};
  ~W() { return; };
  virtual void Add(Int_t m, Int_t n);
  void Set_i(Int_t m) { i = m;};
  void Set_j(Int_t n) { j = n;};
  Int_t i,j;
private:
};
