CREATE TABLE Tests(
    test_id VARCHAR(100) NOT NULL,
    test_min_score INTEGER NOT NULL CHECK (typeof(test_min_score) = 'integer'),
    test_max_score INTEGER NOT NULL CHECK (typeof(test_max_score) = 'integer'),
	  
	PRIMARY KEY(test_id),
	  
	CONSTRAINT valid_test_score_range CHECK(test_max_score > test_min_score)
);--#split#--
	  
	  
CREATE TABLE Clusters(
    test_id VARCHAR(100) NOT NULL,
    cluster_nbr INTEGER NOT NULL CHECK (typeof(cluster_nbr) = 'integer'),
    cluster_name VARCHAR(100) NOT NULL,
      
    PRIMARY KEY  (test_id,cluster_nbr),
	  
	FOREIGN KEY(test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE
   );--#split#--
 
 
 CREATE TABLE Cluster_scores(
    test_id VARCHAR(100) NOT NULL,
    cluster_nbr INTEGER NOT NULL,
    cluster_score INTEGER NOT NULL,
    test_score_est INTEGER NOT NULL, -- estimated test score based on cluster score
      
    PRIMARY KEY (test_id, cluster_nbr, cluster_score),      
    FOREIGN KEY (test_id, cluster_nbr) REFERENCES Clusters(test_id, cluster_nbr) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--

CREATE TABLE Items(
	item_id VARCHAR(50) NOT NULL PRIMARY KEY, -- assessment item identifier in case of qti items, in case of html the file name without extension or path 
	item_name VARCHAR(100) NOT NULL DEFAULT ''
);--#split#--

CREATE TABLE Cluster_design(
	test_id VARCHAR(100) NOT NULL,
	cluster_nbr INTEGER NOT NULL,
	item_nbr INTEGER NOT NULL,
	item_id VARCHAR(50) NOT NULL,	
	
	PRIMARY KEY (test_id,cluster_nbr,item_nbr),
	UNIQUE(test_id,cluster_nbr,item_id),

	FOREIGN KEY	(test_id, cluster_nbr) REFERENCES Clusters(test_id, cluster_nbr) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (item_id) REFERENCES Items (item_id) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--
	

CREATE TABLE Standards(
    test_id VARCHAR(100) NOT NULL,
    standard_nbr INTEGER NOT NULL CHECK (typeof(standard_nbr) = 'integer'),
    standard_name VARCHAR(100) NOT NULL,
      
    PRIMARY KEY (test_id,standard_nbr),
	  
	FOREIGN KEY(test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--


CREATE TABLE Population(
	test_id VARCHAR(100) NOT NULL,
	test_score INTEGER NOT NULL, 	
		/* the total score on the test, it is assumed that this is on the same scale as 
			test_score_est in the Cluster_scores table. */
	test_score_frequency INTEGER NOT NULL,
		-- to save some space and bandwith, this is saved as score plus frequency rather than saving all the individual scores
	
	PRIMARY KEY (test_id, test_score),
 
	FOREIGN KEY (test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE
 );--#split#--
 
 CREATE TABLE Users(
 	username VARCHAR(50) NOT NULL,
	user_password VARCHAR(250) NOT NULL,
	user_role VARCHAR(50) NOT NULL CHECK ( user_role IN('rater','group_leader') ), 
	user_realname VARCHAR(100) NOT NULL DEFAULT '',
	  
	PRIMARY KEY (username),
	UNIQUE (username, user_role)
 
 );--#split#--
 
 CREATE TABLE Group_leaders(
 	  username VARCHAR(50) NOT NULL,
	  user_role VARCHAR(50) NOT NULL DEFAULT 'group_leader' CHECK(user_role = 'group_leader'),	  

	  PRIMARY KEY (username),
	  FOREIGN KEY (username,user_role) REFERENCES Users (username,user_role) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--
 
CREATE TRIGGER insert_group_leader AFTER INSERT ON Users
	WHEN NEW.user_role = 'group_leader'
BEGIN
	INSERT INTO Group_leaders (username) VALUES(NEW.username);
END;--#split#--
 
 CREATE TABLE Group_leader_test_assignment(
	username VARCHAR(50) NOT NULL,
	test_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY(username,test_id),
	FOREIGN KEY(test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY(username) REFERENCES group_leaders(username) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--
	
 
CREATE TABLE Raters(
	username VARCHAR(50) NOT NULL,
	user_role VARCHAR(50) NOT NULL DEFAULT 'rater' CHECK(user_role = 'rater'),	  
	PRIMARY KEY (username),
	FOREIGN KEY(username, user_role) REFERENCES Users(username,user_role) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--


CREATE TRIGGER insert_rater AFTER INSERT ON Users
	WHEN NEW.user_role = 'rater'
BEGIN
	INSERT INTO Raters (username) VALUES(NEW.username);
END;--#split#--

CREATE TABLE Rater_test_assignment(
	username VARCHAR(50) NOT NULL,
	test_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY(username,test_id),
	FOREIGN KEY(test_id) REFERENCES Tests(test_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY(username) REFERENCES Raters(username) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--


CREATE TABLE Standard_setting_design(
	test_id VARCHAR(100) NOT NULL,
    standard_setting_round INTEGER NOT NULL CHECK(standard_setting_round IN(1,2)),
    standard_setting_round_status VARCHAR(10) NOT NULL DEFAULT 'open' CHECK(standard_setting_round_status IN('open','closed')),
    standard_nbr INTEGER NOT NULL,
      
    PRIMARY KEY (test_id, standard_setting_round, standard_nbr),      
    FOREIGN KEY (test_id,standard_nbr) REFERENCES Standards(test_id,standard_nbr) ON UPDATE CASCADE ON DELETE CASCADE
 );--#split#--

CREATE TABLE Standard_setting_scores(
    test_id VARCHAR(100) NOT NULL,
    cluster_nbr INTEGER NOT NULL,
    cluster_score INTEGER NOT NULL,
    standard_nbr INTEGER NOT NULL,
    username VARCHAR(50) NOT NULL,
    standard_setting_round INTEGER NOT NULL,
      
    PRIMARY KEY	(test_id, standard_nbr, cluster_nbr, username, standard_setting_round),
      
    FOREIGN KEY	(test_id, cluster_nbr, cluster_score) REFERENCES Cluster_scores(test_id, cluster_nbr, cluster_score) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY	(test_id, standard_setting_round,standard_nbr) REFERENCES Standard_setting_design(test_id,  standard_setting_round,standard_nbr) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (username,test_id) REFERENCES Rater_test_assignment(username,test_id) ON UPDATE CASCADE ON DELETE CASCADE
);--#split#--

CREATE VIEW Users_test_assignment AS 
	SELECT username, test_id, 'group_leader' AS user_role FROM Group_leader_test_assignment
	UNION
	SELECT username, test_id, 'rater' AS user_role FROM Rater_test_assignment;


