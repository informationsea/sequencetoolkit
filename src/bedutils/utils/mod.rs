use bio::data_structures::interval_tree::*;
use std::collections::HashSet;
use std::ops::Range;

pub struct RegionMerger {
    interval_tree: IntervalTree<u64, (Range<u64>, u32)>,
    region_set: HashSet<Range<u64>>,
}

impl Default for RegionMerger {
    fn default() -> Self {
        RegionMerger::new()
    }
}

impl RegionMerger {
    pub fn new() -> Self {
        RegionMerger {
            interval_tree: IntervalTree::new(),
            region_set: HashSet::new(),
        }
    }
    pub fn add_interval(&mut self, interval: Range<u64>) {
        let mut start: u64 = interval.start;
        let mut end: u64 = interval.end;
        let mut count = 1;
        for one in self.interval_tree.find(interval) {
            let data = one.data();
            start = start.min(data.0.start);
            end = end.max(data.0.end);
            self.region_set.remove(&data.0);
            count += data.1;
        }
        self.interval_tree.insert(start..end, (start..end, count));
        self.region_set.insert(start..end);
    }

    pub fn regions(&self) -> &HashSet<Range<u64>> {
        &self.region_set
    }
}
