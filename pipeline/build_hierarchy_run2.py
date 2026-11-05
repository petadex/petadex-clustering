"""Run 2 hierarchy (2026-11-05): one row per input core, same columns as Run 1's results/hierarchy.tsv
(domain_id petadex_id source component c90 c60 c30) plus cath_domain (Run 2 clusters per CATH domain).
Built from run2/L{90,60,30}/<domain>.all.tsv (rep \t member). Nesting: L60 members are L90 centroids, L30 members
are L60 centroids. Read-only on run2/L*/; writes run2/products/hierarchy_run2.parquet."""
import glob, os, duckdb
R = "/watson/rnalab/dzhu/petadex-clustering/run2"; O = f"{R}/products"
con = duckdb.connect(); con.execute(f"set memory_limit='48GB'; set threads=16; set temp_directory='{O}/tmp'; set preserve_insertion_order=false")
def lvl(L):
    fs = [f for f in sorted(glob.glob(f"{R}/{L}/*.all.tsv"))]
    return " union all ".join(f"""select column0 rep, column1 mem, '{os.path.basename(f)[:-8]}' dom from read_csv('{f}', delim='\t',
        header=false, quote='', escape='', all_varchar=true, columns={{'column0':'VARCHAR','column1':'VARCHAR'}})""" for f in fs)
for L in ("L90", "L60", "L30"):
    if os.path.exists(f"{O}/tmp_{L}.parquet"): continue
    con.execute(f"copy ({lvl(L)}) to '{O}/tmp_{L}.parquet' (format parquet, compression zstd)"); print(L, "read", flush=True)
con.execute(f"""copy (
  select a.mem domain_id, split_part(a.mem, '|', 1) petadex_id,
         case left(a.mem, 2) when 'PL' then 'pazy' when 'NR' then 'nr' when 'LG' then 'logan' end "source",
         split_part(a.mem, '|', 2) component, a.dom cath_domain, a.rep c90, b.rep c60, c.rep c30
  from '{O}/tmp_L90.parquet' a left join '{O}/tmp_L60.parquet' b on b.mem = a.rep
                               left join '{O}/tmp_L30.parquet' c on c.mem = b.rep)
  to '{O}/hierarchy_run2.parquet' (format parquet, compression zstd, row_group_size 1000000)""")
print("written", flush=True)
H = f"{O}/hierarchy_run2.parquet"; q = lambda s: con.execute(s).fetchone()
print("rows", q(f"select count(*), count(distinct domain_id) from '{H}'"))
print("null parents", q(f"""select count_if(c60 is null), count_if(c30 is null), count_if("source" is null) from '{H}'"""))
print("distinct c90/c60/c30", q(f"select count(distinct c90), count(distinct c60), count(distinct c30) from '{H}'"))
print("c90 in one c60 / c60 in one c30", q(f"select count(*) filter (where n>1) from (select c90, count(distinct c60) n from '{H}' group by 1)"),
      q(f"select count(*) filter (where n>1) from (select c60, count(distinct c30) n from '{H}' group by 1)"))
print("by source", con.execute(f"""select "source", count(*), count(distinct petadex_id) from '{H}' group by 1 order by 1""").fetchall())
print("c90 centroid source", con.execute(f"select left(c90,2), count(distinct c90) from '{H}' group by 1 order by 1").fetchall())
for L in ("L90", "L60", "L30"): os.remove(f"{O}/tmp_{L}.parquet")
