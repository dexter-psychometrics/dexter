pragma foreign_keys=1;--#split#--

CREATE TABLE dxItems(
	item_id VARCHAR(100) NOT NULL,
	PRIMARY KEY(item_id)
) WITHOUT ROWID;--#split#--

CREATE TABLE dxScoring_rules(
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	item_score INTEGER NOT NULL, -- might add text, like 'missing', etc.
	
	PRIMARY KEY (item_id, response),
	FOREIGN KEY (item_id) REFERENCES dxItems(item_id) ON UPDATE CASCADE ON DELETE CASCADE
) WITHOUT ROWID;--#split#--


CREATE TABLE dxBooklets(
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id)	
) WITHOUT ROWID;--#split#-- 

CREATE TABLE dxTestparts(
	booklet_id VARCHAR(100) NOT NULL,
	testpart_nbr INTEGER NOT NULL DEFAULT 1 CHECK(testpart_nbr >= 1),
	
	PRIMARY KEY(booklet_id, testpart_nbr)
) WITHOUT ROWID;--#split#--


CREATE TABLE dxBooklet_design(
	booklet_id VARCHAR(100) NOT NULL,
	testpart_nbr INTEGER NOT NULL DEFAULT 1,
	item_id VARCHAR(100) NOT NULL,
	item_position INTEGER NOT NULL CHECK(item_position >= 1),
	
	PRIMARY KEY (booklet_id,testpart_nbr, item_id),
	UNIQUE		(booklet_id,testpart_nbr, item_position),
	
	FOREIGN KEY (item_id) REFERENCES dxItems(item_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id, testpart_nbr) REFERENCES dxtestparts(booklet_id, testpart_nbr) ON UPDATE CASCADE ON DELETE CASCADE
) WITHOUT ROWID;--#split#--


CREATE TABLE dxMulti_stage_rules(
	booklet_id VARCHAR(100) NOT NULL,
	from_testpart_nbr INTEGER NOT NULL,
	to_testpart_nbr INTEGER NOT NULL,
	from_testpart_min_score INTEGER NOT NULL CHECK(from_testpart_min_score >= 0),
	from_testpart_max_score INTEGER NOT NULL CHECK(from_testpart_max_score >= 0),	-- both min and max are inclusive
																					-- max is intentionally allowed to be set to larger than the possible score on a testpart
	
	PRIMARY KEY (booklet_id, from_testpart_nbr, to_testpart_nbr),
	
	FOREIGN KEY (booklet_id, from_testpart_nbr) REFERENCES dxTestparts(booklet_id, testpart_nbr) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id, to_testpart_nbr) REFERENCES dxTestparts(booklet_id, testpart_nbr) ON UPDATE CASCADE ON DELETE CASCADE,

	CONSTRAINT move_forward CHECK(to_testpart_nbr > from_testpart_nbr),	-- this also implies there cannot be a rule to go to the first testpart
	CONSTRAINT max_score_gte_min_score CHECK(from_testpart_max_score >= from_testpart_min_score)
) WITHOUT ROWID;--#split#--	


 
CREATE TABLE dxPersons(
	person_id VARCHAR(100) NOT NULL PRIMARY KEY
)  WITHOUT ROWID;--#split#--



CREATE TABLE dxAdministrations(
	person_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id, person_id),

	FOREIGN KEY (person_id) REFERENCES dxPersons(person_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id) REFERENCES dxBooklets(booklet_id) ON UPDATE CASCADE ON DELETE CASCADE
) WITHOUT ROWID;--#split#--

CREATE TABLE dxResponses(
	person_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	testpart_nbr INTEGER NOT NULL DEFAULT 1,
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id, person_id,  item_id),
	
	FOREIGN KEY (booklet_id, person_id) REFERENCES dxAdministrations(booklet_id, person_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (booklet_id, testpart_nbr, item_id) REFERENCES dxBooklet_design(booklet_id, testpart_nbr, item_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (item_id, response) REFERENCES dxScoring_rules(item_id, response) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED
	-- foreign key constraints deferred for speed
) WITHOUT ROWID;--#split#--

CREATE VIEW dxBooklet_stats AS
WITH 
A1 AS ( SELECT booklet_id, COUNT(*) AS n_persons FROM dxAdministrations GROUP BY booklet_id),
I1 AS ( SELECT booklet_id, COUNT(DISTINCT item_id) AS n_items FROM  dxBooklet_design  GROUP BY booklet_id)
SELECT B.booklet_id, COALESCE(n_items,0) AS n_items, COALESCE(n_persons,0) AS n_persons
	FROM dxBooklets AS B
		LEFT OUTER JOIN A1 USING(booklet_id)
			LEFT OUTER JOIN I1 USING(booklet_id);--#split#--


-- constraint to make sure there can be no overlap, implement check to ensure this query returns no rows.

CREATE TRIGGER ms_rules_rng_overlap_ins AFTER INSERT ON dxMulti_stage_rules
BEGIN
	SELECT RAISE(ROLLBACK, 'the range [min_score, max_score] may not overlap with another rule on the same testpart within the same test')
		WHERE EXISTS(
			SELECT booklet_id, from_testpart_nbr
				FROM dxMulti_stage_rules AS M1
					INNER JOIN  dxMulti_stage_rules AS M2 
														USING(booklet_id, from_testpart_nbr)
				WHERE M1.from_testpart_min_score BETWEEN M2.from_testpart_min_score AND M2.from_testpart_max_score 
						OR M1.from_testpart_max_score BETWEEN M2.from_testpart_min_score AND M2.from_testpart_max_score
				
				GROUP BY booklet_id, from_testpart_nbr
					HAVING COUNT(*) > 1);
END;
--#split#--

CREATE TRIGGER ms_rules_rng_overlap_upd AFTER UPDATE ON dxMulti_stage_rules
BEGIN
	SELECT RAISE(ROLLBACK, 'the range [min_score, max_score] may not overlap with anotrer rule on the same testpart within the same test')
		WHERE EXISTS(
			SELECT booklet_id, from_testpart_nbr
				FROM dxMulti_stage_rules AS M1
					INNER JOIN dxMulti_stage_rules AS M2 
														USING(booklet_id, from_testpart_nbr)
				WHERE M1.from_testpart_min_score BETWEEN M2.from_testpart_min_score AND M2.from_testpart_max_score 
						OR M1.from_testpart_max_score BETWEEN M2.from_testpart_min_score AND M2.from_testpart_max_score
				
				GROUP BY booklet_id, from_testpart_nbr
					HAVING COUNT(*) > 1);
END;
