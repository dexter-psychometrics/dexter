

CREATE TABLE dxItems(
	item_id VARCHAR(100) NOT NULL,
	PRIMARY KEY(item_id)
) ;--#split#--

CREATE TABLE dxScoring_rules(
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	item_score INTEGER NOT NULL, 
	
	PRIMARY KEY (item_id, response),
	FOREIGN KEY (item_id) REFERENCES dxItems(item_id) ON UPDATE CASCADE ON DELETE CASCADE
) ;--#split#--


CREATE TABLE dxBooklets(
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id)	
) ;--#split#-- 

CREATE TABLE dxTestparts(
	booklet_id VARCHAR(100) NOT NULL,
	testpart_nbr INTEGER NOT NULL DEFAULT 1 CHECK(testpart_nbr >= 1),
	
	PRIMARY KEY(booklet_id, testpart_nbr)
) ;--#split#--


CREATE TABLE dxBooklet_design(
	booklet_id VARCHAR(100) NOT NULL,
	testpart_nbr INTEGER NOT NULL DEFAULT 1,
	item_id VARCHAR(100) NOT NULL,
	item_position INTEGER NOT NULL CHECK(item_position >= 1),
	
	PRIMARY KEY (booklet_id, testpart_nbr, item_id),
	UNIQUE		(booklet_id, testpart_nbr, item_position),
	
	FOREIGN KEY (item_id) REFERENCES dxItems(item_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id, testpart_nbr) REFERENCES dxtestparts(booklet_id, testpart_nbr) ON UPDATE CASCADE ON DELETE CASCADE
) ;--#split#--


CREATE TABLE dxMulti_stage_rules(
	booklet_id VARCHAR(100) NOT NULL,
	from_testpart_nbr INTEGER NOT NULL,
	to_testpart_nbr INTEGER NOT NULL,
	from_testpart_min_score INTEGER NOT NULL CHECK(from_testpart_min_score >= 0),
	from_testpart_max_score INTEGER NOT NULL CHECK(from_testpart_max_score >= 0),	
	-- both min and max are inclusive
	-- max is intentionally allowed to be set to larger than the possible score on a testpart
	
	PRIMARY KEY (booklet_id, from_testpart_nbr, to_testpart_nbr),
	
	FOREIGN KEY (booklet_id, from_testpart_nbr) REFERENCES dxTestparts(booklet_id, testpart_nbr) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id, to_testpart_nbr) REFERENCES dxTestparts(booklet_id, testpart_nbr) ON UPDATE CASCADE ON DELETE CASCADE,

	CONSTRAINT move_forward CHECK(to_testpart_nbr > from_testpart_nbr),	
	-- this also implies there cannot be a rule to go to the first testpart
	CONSTRAINT max_score_gte_min_score CHECK(from_testpart_max_score >= from_testpart_min_score)
) ;--#split#--	


 
CREATE TABLE dxPersons(
	person_id VARCHAR(100) NOT NULL PRIMARY KEY
)  ;--#split#--



CREATE TABLE dxAdministrations(
	person_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id, person_id),

	FOREIGN KEY (person_id) REFERENCES dxPersons(person_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id) REFERENCES dxBooklets(booklet_id) ON UPDATE CASCADE ON DELETE CASCADE
) ;--#split#--

CREATE TABLE dxResponses(
	person_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	testpart_nbr INTEGER NOT NULL DEFAULT 1,
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id, person_id,  item_id),
	
	-- foreign key constraints deferred for speed of insertion
	FOREIGN KEY (booklet_id, person_id) REFERENCES dxAdministrations(booklet_id, person_id) 
	      ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (booklet_id, testpart_nbr, item_id) REFERENCES dxBooklet_design(booklet_id, testpart_nbr, item_id) 
	      ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (item_id, response) REFERENCES dxScoring_rules(item_id, response) 
	      ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED
) ;--#split#--

CREATE VIEW dxBooklet_stats AS
WITH 
A1 AS ( SELECT booklet_id, COUNT(*) AS n_persons FROM dxAdministrations GROUP BY booklet_id),
I1 AS ( SELECT booklet_id, COUNT(DISTINCT item_id) AS n_items FROM  dxBooklet_design  GROUP BY booklet_id)
SELECT B.booklet_id, COALESCE(n_items,0) AS n_items, COALESCE(n_persons,0) AS n_persons
	FROM dxBooklets AS B
		LEFT OUTER JOIN A1 USING(booklet_id)
			LEFT OUTER JOIN I1 USING(booklet_id);
			
--triggers on multi stage rules are omitted since practically no enigne supports standard sql for triggers