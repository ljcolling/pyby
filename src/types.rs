use dict_derive::{FromPyObject, IntoPyObject};
use itertools::Itertools;
use itertools::izip;

#[derive(FromPyObject, IntoPyObject, Debug)]
pub struct SamplingRule {
    n_min: i32,
    n_max: i32,
    step_size: i32,
}


impl SamplingRule {
    pub fn unpack(&self) -> (i32, i32, i32) {
        (self.n_min, self.n_max, self.step_size)
    }
}



impl SimTable {
    pub fn table(self) -> Vec<(usize, f64, i32, Option<f64>, f64, f64, f64)> {
        let mut id = Vec::new();
        let mut effsize = Vec::new();
        let mut n = Vec::new();
        let mut bf = Vec::new();
        let mut emp_es = Vec::new();
        let mut t = Vec::new();
        let mut pvalue = Vec::new();

        self.rows.into_iter().for_each(|x| {
            id.push(x.id);
            effsize.push(x.effsize);
            n.push(x.n);
            bf.push(x.bf);
            emp_es.push(x.emp_es);
            t.push(x.t);
            pvalue.push(x.pvalue);
        });

        izip!(id, effsize, n, bf, emp_es, t, pvalue).collect_vec()
    }
}

#[derive(Debug)]
pub struct SimRow {
   pub id: usize,
   pub effsize: f64,
   pub n: i32,
   pub bf: Option<f64>,
   pub emp_es: f64,
   pub t: f64,
   pub pvalue: f64,
}

pub struct SimTable {
    pub rows: Vec<SimRow>,
}

#[derive(Debug)]
struct PairedSim {
    id: Vec<usize>,
    effsize: Vec<f64>,
    n: Vec<i32>,
    bf: Vec<Option<f64>>,
    emp_es: Vec<f64>,
    t: Vec<f64>,
    pvalue: Vec<f64>,
}
impl PairedSim {
    fn new() -> PairedSim {
        PairedSim {
            id: Vec::new(),
            effsize: Vec::new(),
            n: Vec::new(),
            bf: Vec::new(),
            emp_es: Vec::new(),
            t: Vec::new(),
            pvalue: Vec::new(),
        }
    }
}
impl FromIterator<(usize, f64, i32, Option<f64>, f64, f64, f64)> for PairedSim {
    fn from_iter<I: IntoIterator<Item = (usize, f64, i32, Option<f64>, f64, f64, f64)>>(
        iter: I,
    ) -> Self {
        let mut c = PairedSim::new();

        for item in iter {
            c.id.push(item.0);
            c.effsize.push(item.1);
            c.n.push(item.2);
            c.bf.push(item.3);
            c.emp_es.push(item.4);
            c.t.push(item.5);
            c.pvalue.push(item.6);
        }

        c
    }
}

